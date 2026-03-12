-- EquivalenceWorker.m2
--
-- Robust parallel equivalence computation for CY topology classification.
-- Multiple M2 processes can run this simultaneously.
-- Safe to kill at any time: unfinished items stay in inprogress/ for recovery.
--
-- Usage:
--   needsPackage "StringTorics"
--   load "EquivalenceWorker.m2"
--   setupWorkDir "./work"
--
--   -- Create todo items from your data:
--   createTodoItems(yourItemList, "h11-4")  -- prefix for filenames
--
--   -- Load your CY data:
--   Xs = hashTable { ... };  -- label => CalabiYauInToric
--
--   -- Run the worker:
--   runWorker(Xs, findEquivalence)
--
--   -- Or with a timeout (in seconds):
--   runWorker(Xs, findEquivalence, TimeLimit => 120)
--
--   -- Recovery after crash (moves inprogress items back to todo):
--   recoverItems()
--
-- File format:
--   Each todo/done file contains a valid M2 expression: a list of sublists.
--   Each sublist is {repLabel, {label2, matrix2}, {label3, matrix3}, ...}
--   where repLabel is a Sequence (i,j) and matrices are over ZZ.
--   Elements within a sublist are equivalent, witnessed by the matrix
--   (which maps the element to the representative).

-- =============================================
-- Directory setup
-- =============================================
WorkDir = null
TodoDir = null
InProgressDir = null
DoneDir = null

setupWorkDir = (dir) -> (
    WorkDir = dir;
    TodoDir = dir | "/todo";
    InProgressDir = dir | "/inprogress";
    DoneDir = dir | "/done";
    run("mkdir -p '" | TodoDir | "'");
    run("mkdir -p '" | InProgressDir | "'");
    run("mkdir -p '" | DoneDir | "'");
)

-- =============================================
-- File I/O
-- =============================================

-- Write a partition (list of sublists) to a file.
-- Each sublist: {repLabel, {label, matrix}, ...}
writePartition = (filename, partition) -> (
    f := openOut filename;
    f << "{" << endl;
    for i from 0 to #partition - 1 do (
        sub := partition#i;
        f << "  {" << toString sub#0;
        for j from 1 to #sub - 1 do (
            f << "," << endl;
            f << "    {" << toString(sub#j#0) << ", " << toString(sub#j#1) << "}";
        );
        f << "}";
        if i < #partition - 1 then f << ",";
        f << endl;
    );
    f << "}" << endl;
    close f;
)

-- Read a partition from a file.
readPartition = (filename) -> value get filename

-- =============================================
-- Claiming items (atomic via mv)
-- =============================================

-- Try to claim the next available todo item.
-- Returns (filename, filepath) on success, or null if no items remain.
claimItem = () -> (
    files := sort select(readDirectory TodoDir, f -> match("\\.m2$", f));
    for f in files do (
        src := TodoDir | "/" | f;
        dst := InProgressDir | "/" | f;
        if 0 == run("mv '" | src | "' '" | dst | "'") then (
            return (f, dst);
        );
        -- If mv failed, another process claimed it; try next
    );
    null
)

-- =============================================
-- Processing logic
-- =============================================

-- Process a single item file.
-- Returns (refinedPartition, indeterminatePairs).
--
-- findEquiv should be a function (label1, label2, Xs) -> (status, data)
-- where status is CONSISTENT, INCONSISTENT, or INDETERMINATE.
processItem = (filepath, Xs, findEquiv, timeLimit) -> (
    partition := readPartition filepath;
    n := #partition;
    if n <= 1 then return (partition, {});

    -- Mutable copies for merging
    groups := new MutableList from partition;
    alive := new MutableList from toList(n : true);
    indeterminate := {};

    for i from 0 to n - 1 do (
        if not alive#i then continue;
        repI := (groups#i)#0;
        for j from i + 1 to n - 1 do (
            if not alive#j then continue;
            repJ := (groups#j)#0;

            << "    comparing " << toString repI << " vs " << toString repJ << "..." << flush;

            result := null;
            if timeLimit > 0 then (
                try alarm timeLimit;
                try (
                    result = findEquiv(repI, repJ, Xs);
                ) then (
                    alarm 0;  -- cancel alarm
                ) else (
                    alarm 0;
                    << " TIMEOUT/ERROR" << endl;
                    indeterminate = append(indeterminate, {repI, repJ});
                    continue;
                );
            ) else (
                try (
                    result = findEquiv(repI, repJ, Xs);
                ) else (
                    << " ERROR" << endl;
                    indeterminate = append(indeterminate, {repI, repJ});
                    continue;
                );
            );

            status := result#0;
            if status === CONSISTENT then (
                << " CONSISTENT" << endl;
                M := result#1;
                -- Merge group j into group i.
                -- Each element {label, Mk} in group j satisfies: label ~ repJ via Mk.
                -- We have repI ~ repJ via M.
                -- So label ~ repI via M * Mk.
                gj := groups#j;
                -- Add repJ itself with matrix M
                newEntries := {{repJ, transpose M}};
                for k from 1 to #gj - 1 do (
                    newEntries = append(newEntries, {gj#k#0, gj#k#1 * (transpose M^-1)});
                );
                groups#i = join(groups#i, newEntries);
                alive#j = false;
            )
            else if status === INDETERMINATE then (
                << " INDETERMINATE" << endl;
                indeterminate = append(indeterminate, {repI, repJ});
            )
            else (
                -- INCONSISTENT
                << " INCONSISTENT" << endl;
            );
        );
    );

    -- Collect surviving groups
    refinedPartition := for i from 0 to n - 1 list (
        if alive#i then groups#i else continue
    );
    (refinedPartition, indeterminate)
)

-- =============================================
-- Worker loop
-- =============================================

runWorker = {TimeLimit => 0} >> opts -> (Xs, findEquiv) -> (
    if WorkDir === null then error "Call setupWorkDir first";
    count := 0;
    timeLimit := opts.TimeLimit;
    while true do (
        claimed := claimItem();
        if claimed === null then (
            << "No more items. Processed " << count << " items." << endl;
            break;
        );
        (itemName, itemPath) := claimed;
        << "Processing " << itemName << "..." << endl;

        (result, indeterminate) := processItem(itemPath, Xs, findEquiv, timeLimit);

        -- Write result to done/
        donePath := DoneDir | "/" | itemName;
        writePartition(donePath, result);

        -- Write indeterminate pairs if any
        if #indeterminate > 0 then (
            indPath := DoneDir | "/" | replace("\\.m2$", "-indeterminate.txt", itemName);
            f := openOut indPath;
            for p in indeterminate do (
                f << toString(p#0) << "  " << toString(p#1) << endl;
            );
            close f;
        );

        -- Remove from inprogress (item is safely in done/ now)
        removeFile itemPath;

        nGroups := #result;
        << "  Done: " << itemName << " -> " << nGroups << " groups";
        if #indeterminate > 0 then << " (" << #indeterminate << " indeterminate)";
        << endl;
        count = count + 1;
    );
)

-- =============================================
-- Recovery and utilities
-- =============================================

-- Move all inprogress items back to todo (run after a crash).
recoverItems = () -> (
    files := select(readDirectory InProgressDir, f -> match("\\.m2$", f));
    for f in files do (
        run("mv '" | InProgressDir | "/" | f | "' '" | TodoDir | "/" | f | "'");
    );
    << "Recovered " << #files << " items." << endl;
)

-- Create todo items from a list of partitions.
-- items: a list of pairs (name, partition), where name is a string/number
--        and partition is a list of sublists in the format above.
-- prefix: string prefix for filenames.
createTodoItems = (items, prefix) -> (
    if WorkDir === null then error "Call setupWorkDir first";
    for item in items do (
        (name, partition) := item;
        filename := TodoDir | "/" | prefix | "-" | toString name | ".m2";
        writePartition(filename, partition);
    );
    << "Created " << #items << " todo items." << endl;
)

-- Show status of the work queue.
workStatus = () -> (
    if WorkDir === null then error "Call setupWorkDir first";
    todoFiles := select(readDirectory TodoDir, f -> match("\\.m2$", f));
    inprogFiles := select(readDirectory InProgressDir, f -> match("\\.m2$", f));
    doneFiles := select(readDirectory DoneDir, f -> match("\\.m2$", f));
    << "Todo:       " << #todoFiles << endl;
    << "InProgress: " << #inprogFiles << endl;
    << "Done:       " << #doneFiles << endl;
)

-- Collect all results from done/ into a single list.
collectResults = () -> (
    if WorkDir === null then error "Call setupWorkDir first";
    files := sort select(readDirectory DoneDir, f -> match("\\.m2$", f) and not match("indeterminate", f));
    for f in files list (
        (f, readPartition(DoneDir | "/" | f))
    )
)

end--

restart
needs "EquivalenceWorker.m2"
setupWorkDir "./work-h11-5"
