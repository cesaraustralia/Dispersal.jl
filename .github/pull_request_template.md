# Description

[Describe what this PR adds to the package]

## Changes (choose one only, with x or click after creating PR):

- [ ] New formulation/rule
- [ ] Documentation
- [ ] Bug fix for existing rules (non-breaking syntax)
- [ ] Minor change to existing rules (no breaks to old syntax)
- [ ] Major change to existing rules (breaking syntax, needs a breaking version bump)

# Checklist:

- [ ] My code follows the [BlueStyle](https://github.com/invenia/BlueStyle) style guide
- [ ] My code is commented with explanations of anything hard to understand
- [ ] There are no unnecessary or unrelated code changes included in this PR
- [ ] There are no dependencies added: if so, please consider if they are necessary and explain why they are required

### For bugfixes

- [ ] My PR addresses a single bug/connected group of bugs
- [ ] I have included a test that will ensure the bug does not reoccur

### For new `Rule`s

- [ ] My PR is limited to a clear theme/conceptual unit (otherwise make multiple PRs that can be discussed separately).
- [ ] New rules have tests in a file of the same name as the rule file, e.g. `/test/newrulename.jl`
- [ ] I have cited the source of the new formulation/s, if applicable
- [ ] I have added a documentation string to formulation struct/s, and added an entry in the docs
- [ ] The doc string has an Example block using `jldoctest`, like :
    ````
    """
        TheRule <: CellRule
    
    Rest of the docs...
    
    # Example
    
    In this example we do...
    ```jldoctest
    rule = TheRule(; params...)
    output = SomeOutput(; init=the_init_grid, tspan=the_tspan)
    sim!(output, rule);
    ```
    """
    ````
