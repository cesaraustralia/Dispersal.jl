# Description

[Describe what this PR adds to the package]

## Changes (choose one only):

- [ ] New formulation/rule
- [ ] Bug fix for existing rules (non-breaking syntax)
- [ ] Major change to existing rules (breaking syntax, needs a breaking version bump)

# Checklist:

- [ ] My code follows the [BlueStyle](https://github.com/invenia/BlueStyle) style guide
- [ ] Code is commented with explanations of anything hard to understand
- [ ] There are no unnecessary or unrelated code changes included in this PR
- [ ] There are no depencecies added, if so please consider if they are necessary and explain why they are required

### For bugfixes

- [ ] My PR addresses a single bug/connected group of bugs
- [ ] I have included a test that will ensure the bug does not recur

### For new `Rule`s

- [ ] My PR is limited to a clear theme/cenceptual unit. (Otherwise make multiple PRs that can be discussed separately)
- [ ] New rules have tests in a file of the same name as the rule file, e.g. `/test/newrulename.jl`
- [ ] I have cited the source of the new formulation/s, if applicable
- [ ] I have added a documentation string to formulation struct/s, and added an entry in the docs

