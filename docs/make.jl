using Documenter
using NCutSegmentation

makedocs(
    sitename = "NCutSegmentation",
    format = Documenter.HTML(),
    modules = [NCutSegmentation]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
