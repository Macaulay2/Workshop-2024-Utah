-- Scratch code from ReflexivePolytopes.m2
-- Old generateTriangulations implementation via bistellar flips (marked FIX)

-*
generateTriangulations Triangulation := opts -> T -> (
    allT := new MutableHashTable;
    allT#T = true;
    TODO := {T};
    while #TODO > 0 and #(keys allT) < opts.Limit do (
        nextTRI := TODO#0;
        TODO = drop(TODO,1);
        --flips := affineCircuits nextTRI;
        --flips := select(affineCircuits nextTRI, z -> #z#0 > 1 and #z#1 > 1);
        flips := select(affineCircuits nextTRI, z -> #z#0 > 1 or #z#1 > 1);
        fliptris := for f in flips list bistellarFlip(nextTRI, f);
        newT := select(fliptris, x -> x =!= null);
        for T in newT do (
            if not allT#?T then (
                --<< "new triangulation: " << T << endl;
                if not opts.RegularOnly or isRegularTriangulation T then (
                    allT#T = true;
                    TODO = append(TODO, T);
                    );
                ));
        if debugLevel > 0 then
            << "todo = " << #TODO << " and #triang = " << #(keys allT) << endl;
        );
    keys allT
    )

generateTriangulations(Matrix, List) := opts -> (Amat, tri) -> (
    allT := new MutableHashTable;
    allT#tri = true;
    TODO := {tri};
    while #TODO > 0 and #(keys allT) < opts.Limit do (
        nextTRI := TODO#0;
        TODO = drop(TODO,1);
        --flips := affineCircuits nextTRI;
        --flips := select(affineCircuits nextTRI, z -> #z#0 > 1 and #z#1 > 1);
        flips := select(affineCircuits(Amat, tri), z -> #z#0 > 1 and #z#1 > 1);
        fliptris := for f in flips list bistellarFlip(nextTRI, f);
        newT := select(fliptris, x -> x =!= null);
        for T in newT do (
            if not allT#?T then (
                --<< "new triangulation: " << T << endl;
                if not opts.RegularOnly or isRegularTriangulation(Amat, T, Homogenize => false) then (
                    allT#T = true;
                    TODO = append(TODO, T);
                    );
                ));
        if debugLevel > 0 then
            << "todo = " << #TODO << " and #triang = " << #(keys allT) << endl;
        );
    keys allT
    )
*-
