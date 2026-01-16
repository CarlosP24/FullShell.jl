# Package Registration Checklist

Use this checklist to prepare FullShell.jl for registration in the Julia General Registry.

## Pre-Registration Requirements

### Package Structure
- [x] Valid `Project.toml` with name, uuid, version, authors
- [x] Package follows Julia naming conventions (ends in `.jl`)
- [x] Source code in `src/` directory
- [x] Main module file matches package name (`FullShell.jl`)
- [x] Tests in `test/` directory with `runtests.jl`

### Documentation
- [x] All exported functions have docstrings
- [x] README.md with basic usage examples
- [x] LICENSE file
- [x] Documenter.jl setup for hosted documentation
- [ ] Documentation successfully builds and deploys to GitHub Pages
- [x] CONTRIBUTING.md with development guidelines

### Code Quality
- [ ] All tests pass (`Pkg.test("FullShell")`)
- [ ] Package loads without errors
- [ ] Minimal dependencies (review `[deps]` in Project.toml)
- [ ] Compatibility bounds specified for all dependencies
- [x] Julia version compatibility specified (>= 1.9)

### Version Control
- [x] Code hosted on GitHub
- [ ] Clean git history
- [ ] No sensitive information in repository
- [x] Appropriate `.gitignore` file

### Metadata
- [x] Package has a unique UUID
- [ ] Version number follows semantic versioning (currently "1.0.0-DEV")
- [x] Authors listed in Project.toml
- [ ] Package has a DOI (currently: 10.5281/zenodo.11450677)

## Registration Steps

### 1. Finalize Documentation
```bash
# Build documentation locally to verify
cd docs
julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path="..")); Pkg.instantiate()'
julia --project=. make.jl
```

- [x] Documentation builds without errors
- [ ] All API functions documented
- [ ] Examples work as expected
- [ ] Set up GitHub Pages deployment (see DOCUMENTATION.md)

### 2. Set Version Number
In `Project.toml`, change:
```toml
version = "1.0.0-DEV"
```
to:
```toml
version = "0.1.0"  # or "1.0.0" if ready for stable
```

**Version Guidelines:**
- Start with `0.1.0` for initial release (allows API changes)
- Use `1.0.0` only if API is considered stable
- Follow [Semantic Versioning](https://semver.org/):
  - MAJOR: Breaking changes
  - MINOR: New features (backwards compatible)
  - PATCH: Bug fixes

### 3. Run Final Tests
```julia
using Pkg
Pkg.test("FullShell")
```
- [ ] All tests pass
- [ ] No warnings
- [ ] Test coverage is reasonable

### 4. Update README
- [ ] Remove "not yet registered" notice
- [ ] Update installation instructions
- [ ] Verify all links work
- [ ] Add badges (CI, docs, coverage if available)

### 5. Create Git Tag
```bash
git tag -a v0.1.0 -m "Initial release"
git push origin v0.1.0
```

### 6. Register Package
Use the Julia Registrator:

**Option A: Via GitHub Comment**
1. Create a GitHub release for your version tag
2. Comment on a commit/PR: `@JuliaRegistrator register`
3. Follow the bot's instructions

**Option B: Via Web Interface**
1. Go to https://juliahub.com/ui/Registrator
2. Follow the registration wizard

**Option C: Via Julia REPL**
```julia
using LocalRegistry
register("FullShell")  # If you have a local registry setup
```

### 7. Wait for Review
- Registration bot will open a PR to General registry
- Automated checks will run
- May need to make adjustments based on feedback
- Usually takes 1-3 days for approval

## Post-Registration

### After Merge
- [ ] Verify package can be installed: `Pkg.add("FullShell")`
- [ ] Update README with final installation instructions
- [ ] Announce release (if desired)
- [ ] Consider registering on JuliaHub

### Maintenance
- [ ] Set up CI/CD (GitHub Actions) for testing
- [ ] Consider adding code coverage reporting
- [ ] Plan for future releases
- [ ] Monitor issues and PRs

## Common Registration Issues

### Failed Checks
- **Version conflict**: Ensure version number isn't already used
- **Compatibility issues**: Add compat entries for all dependencies
- **Name conflict**: Package name must be unique
- **UUID issues**: UUID must be unique and not change

### How to Fix
1. Make corrections in your repository
2. Update the version number
3. Create a new tag
4. Try registration again

## Additional Resources

- [Pkg.jl Documentation](https://pkgdocs.julialang.org/)
- [General Registry](https://github.com/JuliaRegistries/General)
- [Registrator.jl](https://github.com/JuliaRegistries/Registrator.jl)
- [Julia Package Naming Guidelines](https://pkgdocs.julialang.org/v1/creating-packages/#Package-naming-guidelines)

## Questions?

If you encounter issues during registration:
1. Check the [Pkg.jl documentation](https://pkgdocs.julialang.org/)
2. Ask on [Julia Discourse](https://discourse.julialang.org/)
3. Look at recently registered packages for examples
4. Review the General registry guidelines

---

**Note**: This checklist is based on current Julia package registry requirements as of 2024. 
Requirements may change; always check the official documentation before registering.
