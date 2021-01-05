using DemoActions
using Documenter

makedocs(;
    modules=[DemoActions],
    authors="anand <anandj@uchicago.edu> and contributors",
    repo="https://github.com/anandijain/DemoActions.jl/blob/{commit}{path}#L{line}",
    sitename="DemoActions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://anandijain.github.io/DemoActions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/anandijain/DemoActions.jl",
)
