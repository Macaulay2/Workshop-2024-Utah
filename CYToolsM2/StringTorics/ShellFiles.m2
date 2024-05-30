blocks = method()
blocks(ZZ, ZZ, ZZ) := (ntotal, nperrun, nperfile) -> (
    -- returns list of (filenum, lo, hi), -- starting at lo=0, hi=ntotal-1.
    count := 0;
    filenum := 0;
    thisfile := nperfile;
    while count < ntotal list (
	if thisfile == nperfile then (
	    filenum = filenum + 1;
	    thisfile = 0;
	    );
	lo := count;
	hi := min(ntotal-1, lo + nperrun-1);
	count = hi+1;
	thisfile = thisfile + 1;
	{filenum, lo, hi}
	)
    )

blocks(ZZ, ZZ, ZZ, ZZ) := (abslo, hi, nperrun, nperfile) -> (
    -- returns list of (filenum, lo, hi), -- starting at lo=0, hi=ntotal-1.
    count := 0;
    filenum := 0;
    thisfile := nperfile;
    ntotal := hi-abslo+1;
    while count < ntotal list (
	if thisfile == nperfile then (
	    filenum = filenum + 1;
	    thisfile = 0;
	    );
	lo := count;
	hi = min(ntotal-1, lo + nperrun-1);
	count = hi+1;
	thisfile = thisfile + 1;
	{filenum, lo + abslo, hi + abslo}
	)
    )

createShellFiles = method(Options => {
	Processes => null, -- number of different M2 calls in each file
	Count => null,  -- number of elements per M2 run.
	Range => {0, -1}
	})
createShellFiles(String, String) := opts -> (prefixstr, m2str) -> (
    nperrun := opts.Count;
    nperfile := opts.Processes;
    lo := opts.Range#0;
    hi := opts.Range#1;
    groups := blocks(lo, hi, nperrun, nperfile);
    -- each file will have nperfile lines, plus an echo at the end?
    thisfile := -1;
    F := null;
    for b in groups do (
	if b#0 != thisfile then (
	    if F =!= null then close F;
	    F = openOut(prefixstr|"shell-set-"|b#0|".sh");
	    thisfile = b#0;
	    );
	lo1 := b#1;
	hi1 := b#2;
	str := replace("@lo@", toString lo1, m2str);
	    str =  replace("@hi@", toString hi1, str);
	    str = replace("@prefix@", toString prefixstr, str);
	    F << str << " &" << endl;
	);
    if F =!= null and isOpen F then close F;
    )


end--
  restart
  needs "../ShellFiles.m2"
  m2str = ///M2 -e 'needsPackage "StringTorics"' -e 'createCYDatabase("@prefix@",3,{@lo@,@hi@})' -e 'exit 0' ///
  m2str = ///M2 -e 'needs "ShellFiles.m2"' -e 'testit(@prefix@,{@lo@,@hi@})' -e 'exit 0' ///

blocks(224, 10, 12)  
blocks(0, 223, 10, 12) 
netList blocks(0, 37123, 1000, 10)
createShellFiles("./Foo3/", "callme(@prefix@, @lo@, @hi@)",
    Processes => 10, Count => 1000, Range => {0,17000})

  m2str = ///M2 --stop -e 'needsPackage "StringTorics"' -e 'createCYDatabase("@prefix@",3,{@lo@,@hi@})' -e 'exit 0' ///
createShellFiles("./Foo3/", m2str,
    Processes => 10, Count => 5, Range => {0,243})
