ToricLinearSeries = new Type of HashTable
ToricLinearSeries.synonym = "toric linear series"

toricLinearSeries = method(TypicalValue => ToricLinearSeries)
toricLinearSeries List := ToricLinearSeries => m -> (
    s := new ToricLinearSeries from {
        "monomials" => m  
    };
    s
)

monomials ToricLinearSeries := List => o -> s -> (
    s#"monomials"
)

