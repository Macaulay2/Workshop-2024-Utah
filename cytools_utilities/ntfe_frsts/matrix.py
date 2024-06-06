# -----------------------------------------------------------------------------
# Description:  Simple utility functions for matrix operations.
# -----------------------------------------------------------------------------

import numpy as np
import copy

# basics
# ------
def to_base10(c, B):
    """
    **Description:**
    Converts a number given in components w.r.t. some bases to a number in base
    10.

    **Arguments:**
    - `c` *(array_like)*: A list of the components.
    - `B` *(array_like)*: A list of the bases.

    **Returns:**
    *(numeric)*: The number in base-10
    """
    result = 0
    multiplier = 1
    for c_i, B_i in zip(reversed(c), reversed(B)):
        result += int(c_i) * multiplier
        multiplier *= B_i
    return result

def from_base10(n, B):
    """
    **Description:**
    Split a number in base 10 to components components w.r.t. some bases.

    **Arguments:**
    - `n` *(numeric)*: The number in base 10.
    - `B` *(array_like)*: A list of the bases.

    **Returns:**
    *(list)*: The bases
    """
    c = []
    for B_i in reversed(B):
        c.append(n % B_i)
        n //= B_i
    return list(reversed(c))

# classes
# =======
class LIL():
    """
    This class describes a 2D LIL matrix. This has the same/less functionality
    as scipy.sparse.lil_array, but it is sometimes (much) quicker

    **Arguments:**
    - `width` *(int, optional)*: The width of the matrix.
    - `nworkers` *(integer, optional)*: Number of processors to use for
        multithreaded calculations.
    """
    def __init__(self, dtype, width=None, nworkers=1, iter_densely=False):
        self.arr = []
        self.dtype = dtype
        self.arr_dense = None
        self._len = None
        self.width = width
        self.default_val = 0

        self.iter_densely = iter_densely
        self.nworkers = nworkers

        self._sum_all = None
        self._sum_0 = None
        self._sum_0_dense = None
        self._sum_1 = None

    # basic interface
    # ---------------
    def __repr__(self):
        # piggy-back printing from list
        return self.arr.__repr__()
    
    def __str__(self):
        # piggy-back string conversion from list
        return self.arr.__str__()

    def __iter__(self):
        # iterator
        if self.iter_densely:
            return iter(self.dense())
        else:
            return iter(self.arr)

    def __setitem__(self, idx, value):
        # item assignment
        if not isinstance(idx,tuple):
            raise ValueError(f"Index must be tuple but was {type(idx)}...")

        self.arr[idx[0]][idx[1]] = value

    def __getitem__(self, idx):
        # indexing
        if isinstance(idx,tuple):
            # get element self.arr[i][j]
            if self.width is None:
                print("LIL: Width not set. Inferring from non-zero values...")
                self.width = self.infer_width()

            if idx[1]>=self.width:
                raise IndexError("list index out of range")
            else:
                return self.arr[idx[0]].get(idx[1],0)
        else:
            # get element self.arr[i]
            return self.arr[idx]

    def __len__(self):
        # length
        if self._len is None:
            self._len = len(self.arr)
        return self._len

    def __array__(self, dtype=None):
        # np.array
        return np.array(self.dense(), dtype=dtype)

    @property
    def shape(self):
        return (len(self),self.width)

    def __add__(self, other):
        # addition
        out = LIL(dtype=self.dtype, width=self.width, nworkers=self.nworkers)
        out.arr = self.arr+other.arr
        return out

    # basic methods
    # --------------
    def infer_width(self):
        """
        **Description:**
        Find the minimum width necessary to hold array

        **Arguments:**
        None.
        
        **Returns:**
        Nothing
        """
        return 1+max([max(row.keys()) for row in self.arr])

    def new_row(self):
        """
        **Description:**
        Append an empty row to the dict.

        **Arguments:**
        None.
        
        **Returns:**
        Nothing
        """
        self.arr.append(dict())

    def append(self, toadd, tocopy=True):
        """
        **Description:**
        Append (a) row(s) to the array.

        **Arguments:**
        - `toadd` *(dict or LIL-like)*: Row(s) to add.
        - `tocopy` *(bool,optional)*: Whether to append a copy of toadd.
        
        **Returns:**
        *(LIL)* Itself.
        """
        if len(toadd)==0:
            return self

        # convert to list of dicts
        if isinstance(toadd,dict):
            toadd = [toadd]
        elif isinstance(toadd[0],type(self)):
            toadd = flatten_top(toadd)

        if tocopy:
            self.arr += copy.copy(toadd)
        else:
            self.arr += toadd

        # reset length
        self._len = None
        
        return self

    def col_inds(self):
        return set().union(*[r.keys() for r in self.arr])

    def reindex(self, f=None):
        """
        **Description:**
        Reindex the ith column to be the f(i)-th one.

        **Arguments:**
        - `f` *(dict)*: Dictionary mapping old column indices to new ones.
        
        **Returns:**
        Nothing
        """
        self.arr_dense = None

        # default map is contiguous from 0 to N_cols-1
        if f is None:
            f = {v:k for (k,v) in enumerate(self.col_inds())}

        if self.nworkers<=1:
            for i,row in enumerate(self.arr):
                self.arr[i] = {(f[j] if j in f else j):v for j,v in row.items()}
        else:
            raise NotImplementedError("Multithreading not yet implemented...")

    def unique_rows(self, allow_shuffle=False):
        """
        **Description:**
        Delete repeated rows. Maybe re-orders rows...

        **Arguments:**
        None.
        
        **Returns:**
        Nothing
        """
        self.arr_dense = None

        if allow_shuffle:
            self.arr = [dict(t) for t in {tuple(d.items()) for d in self.arr}]
        else:
            raise NotImplementedError("LIL.unique_rows: Not yet implemented...")

    def dense(self, tocopy=False):
        """
        **Description:**
        Return a dense version of the array

        **Arguments:**
        - `copy` *(bool,optional)*: Whether to return a copy of self.arr_dense.
        
        **Returns:**
        *(np.array)* The dense array
        """
        if self.arr_dense is None:
            # build empty dense array
            height = len(self.arr)

            if self.default_val==0:
                self.arr_dense = np.zeros((height,self.width),dtype=self.dtype)
            else:
                self.arr_dense = self.default_val*np.ones((height,self.width),\
                                                            dtype=self.dtype)

            # fill in output
            if self.nworkers<=1:
                for i,row in enumerate(self.arr):
                    for j,v in row.items():
                        self.arr_dense[i,j] = v
            else:
                raise NotImplementedError("Multithreading not yet "\
                                            "implemented...")

        # return
        if tocopy:
            return self.arr_dense.copy()
        else:
            return self.arr_dense

    def tolist(self):
        return self.dense().tolist()

    def sum(self, axis=None, dense=True):
        if axis is None:
            if self._sum_all is None:
                self._sum_all = np.sum(self.sum(axis=1))
            return self._sum_all
        elif axis==1:
            if self._sum_1 is None:
                self._sum_1 = np.asarray([sum(r.values()) for r in self.arr])
            return self._sum_1
        elif axis==0:
            if dense:
                if self._sum_0_dense is None:
                    self._sum_0_dense = np.asarray([sum(r.get(i,0) for r in\
                                        self.arr) for i in range(self.width)])
                return self._sum_0_dense
            else:
                if self._sum_0 is None:
                    self._sum_0 = {i: sum(r.get(i,0) for r in self.arr) for i\
                                                            in self.col_inds()}
                return self._sum_0

class lazy_tuple:
    """
    A tuple class whose components are only lazily calculated

    **Arguments:**
    - `data` *(misc)*: Tuple elements (or functions to calculate them).
    """
    def __init__(self, *data):
        self._data = tuple(data)
    
    def __repr__(self):
        # piggy-back printing from tuple
        return self._data.__repr__()
    
    def __str__(self):
        # piggy-back string conversion from tuple
        return self._data.__str__()

    def __len__(self):
        # length
        return len(self._data)

    def __getitem__(self, key):
        item = self._data[key]
        
        if (item is None) or callable(item):
            self._data = list(self._data)
            self._data[key] = item()
            self._data = tuple(self._data)
            item = self._data[key]
            
        return item

class LIL_stack():
    """
    This class describes a stack of LIL objects. One could just manually stack the
    rows but this implementation is quicker.

    The stack is organized as a list of options,
        options = [ [top_block_option1, top_block_option2, ...],
                    [next_block_option1, next_block_option2, ...],
                    ...
                    [bot_block_option1, bot_block_option2, ...]]
    and a list of choices
        choices = [i_top_block, i_next_block, ..., i_bot_block]
    E.g., if choices were [7,2,...,6], the stack would look like:
        stack = [top_block_option7;
                 next_block_option2;
                 ...
                 bot_block_option6]

    **Arguments:**
    - `options`: The possible arrays to stack. The entry options[i] is a list
        of all possible blocks you can put in the ith entry
    - `choices`: The selection of which blocks (from options) to stack.
    - `choice_bounds`: The number of possible choices for each block. I.e.,
        [len(opts) for opts in options]
    - `iter_densely`: Whether to iterate densely over the array or sparsely.
    """
    def __init__(self,
                 options: [["ArrayLike"]],
                 choices: [int],
                 choice_bounds: [int],
                 iter_densely: bool = False):

        self._options = options
        if isinstance(choices,int):
            self._choices = choices
        else:
            self._choices = to_base10(choices,choice_bounds)
        self._choice_bounds = choice_bounds
        self.iter_densely = iter_densely

    # basic interfaces
    def __repr__(self):
        # piggy-back printing from list
        return self.arr.__repr__()
    
    def __str__(self):
        # piggy-back string conversion from list
        return self.arr.__str__()

    def __getitem__(self, idx):
        # indexing

        if isinstance(idx,tuple):
            # get element self.arr[i][j]

            if idx[0]>=len(self):
                raise IndexError("0th list index out of range")
            elif idx[0]<0:
                raise IndexError("negative indexing not currently allowed")

            for block in self._blocks():
                L = len(block)
                if idx[0]<L:
                    return block[idx]
                else:
                    idx = (idx[0]-L,idx[1])
        else:
            # get element self.arr[i]
            if idx>=len(self):
                raise IndexError("list index out of range")
            elif idx<0:
                raise IndexError("negative indexing not currently allowed")

            for block in self._blocks():
                L = len(block)
                if idx<L:
                    return block[idx]
                else:
                    idx -= L

    def __len__(self):
        # length
        if not hasattr(self, '_len'):
            self._len = sum(len(opts[i]) for i,opts in
                                            zip(self.choices,self._options))
        return self._len

    def __iter__(self):
        # iterator
        return self._rows(self.iter_densely)

    def __array__(self, dtype=None):
        # np.array
        return np.array(self.dense(), dtype=dtype)

    # properties
    @property
    def choices(self):
        return from_base10(self._choices, self._choice_bounds)
    
    @property
    def dtype(self):
        return self._options[0][0].dtype
    
    @property
    def width(self):
        return self._options[0][0].width

    @property
    def shape(self):
        if not hasattr(self, '_shape'):
            #self._shape = (len(self),self.width) # slow
            self._shape = lazy_tuple(self.__len__,self.width)
        return self._shape
    
    @property
    def is_empty(self):
        if not hasattr(self, '_is_empty'):
            self._is_empty = True
            
            for block in self._blocks():
                if len(block)>0:
                    self._is_empty = False
                    break

        return self._is_empty
    
    def _blocks(self):
        for i,opts in zip(self.choices,self._options):
            yield opts[i]

    def _rows(self,dense=True):
        if dense:
            row_iter = lambda r:r.dense()
        else:
            row_iter = lambda r:r

        for block in self._blocks():
            for h in row_iter(block):
                yield h

    # getter
    @property
    def arr(self):
        if not hasattr(self, '_arr'):
            self._arr = [row for row in self._rows(False)]

        return self._arr

    def dense(self, tocopy=False):
        """
        **Description:**
        Return a dense version of the array

        **Arguments:**
        - `copy` *(bool,optional)*: Whether to return a copy of self.arr_dense.
        
        **Returns:**
        *(np.array)* The dense array
        """
        if not hasattr(self, '_arr_dense'):
            # build empty dense array
            self._arr_dense = np.zeros(self.shape, dtype=self.dtype)

            # fill in output
            for i,row in enumerate(self.arr):
                for j,v in row.items():
                    self._arr_dense[i,j] = v

        # return
        if tocopy:
            return self._arr_dense.copy()
        else:
            return self._arr_dense

    def tolist(self):
        return self.dense().tolist()

    # basic methods
    def sum(self, axis=None, dense=True):
        if axis is None:
            if not hasattr(self, '_sum_all'):
                self._sum_all = np.sum(self.sum(axis=1))
            return self._sum_all
        elif axis==1:
            if not hasattr(self, '_sum_1'):
                self._sum_1 = flatten_top(\
                                    [M.sum(axis=1) for M in self._blocks()],\
                                    as_list=False)
            return self._sum_1
        elif axis==0:
            if dense:
                if not hasattr(self, '_sum_0_dense'):
                    self._sum_0_dense = np.sum(M.sum(axis=0,dense=True) for M\
                                                            in self._blocks())
                return self._sum_0_dense
            else:
                raise NotImplementedError("sparse sum not yet implemented (shouldn't be hard though...)")