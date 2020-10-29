using Weave
# Example on Gradient of Selection
weave("notebook/exampleGradientSelection.jmd";
       doctype = "pandoc2html",
       out_path = "notebook/",
       pandoc_options=["--toc", "-N"])

weave("notebook/exampleGradientSelection.md";
       doctype = "pandoc2html",
       out_path = "notebook/",
       pandoc_options=["--toc", "-N", "--toc-depth=2"])

# Example on Phenology and Dispersal
weave("notebook/examplePhenologyDispersal.jmd";
       doctype = "pandoc2html",
       out_path ="notebook/",
       pandoc_options=["--toc", "-N"])

weave("notebook/examplePhenologyDispersal.md";
       doctype = "pandoc2html",
       out_path = "notebook/",
       pandoc_options=["--toc", "-N", "--toc-depth=2"])