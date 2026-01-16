# Documentation Setup Guide

This guide explains how to use the documentation system that has been set up for FullShell.jl.

## What Has Been Added

### 1. Comprehensive Docstrings
All exported functions and types now have detailed docstrings including:
- Function/type descriptions
- Arguments with types and descriptions
- Return value descriptions
- Usage examples
- Physical interpretations where relevant

### 2. Documentation Structure
A complete Documenter.jl setup has been created:

```
docs/
├── .gitignore           # Ignore build artifacts
├── Project.toml         # Documentation dependencies
├── make.jl              # Documentation build script
└── src/
    ├── index.md         # Main landing page
    ├── examples.md      # Usage examples
    └── api.md           # Complete API reference
```

### 3. GitHub Actions Workflow
Automatic documentation building and deployment via `.github/workflows/documentation.yml`

### 4. Contributing Guidelines
`CONTRIBUTING.md` with instructions for developers

## Building Documentation Locally

### First Time Setup
```bash
cd docs
julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path="..")); Pkg.instantiate()'
```

### Build the Docs
```bash
julia --project=docs docs/make.jl
```

The HTML documentation will be generated in `docs/build/`. Open `docs/build/index.html` in your browser.

### Quick Preview
```bash
cd docs/build
python -m http.server 8000
# Open http://localhost:8000 in your browser
```

## Setting Up GitHub Pages

To enable automatic documentation deployment:

1. **Generate a Documenter key** (one-time setup):
   ```julia
   using DocumenterTools
   DocumenterTools.genkeys(user="CarlosP24", repo="FullShell.jl")
   ```
   This will output a public and private key.

2. **Add the public key to GitHub**:
   - Go to your repository settings → Deploy keys
   - Add a new deploy key with the title "documenter"
   - Paste the public key
   - Check "Allow write access"

3. **Add the private key as a secret**:
   - Go to repository settings → Secrets and variables → Actions
   - Create a new secret named `DOCUMENTER_KEY`
   - Paste the private key

4. **Push your changes** (if not already done):
   ```bash
   git add .
   git commit -m "Add comprehensive documentation"
   git push
   ```

5. **Wait for the GitHub Actions workflow to run**:
   - Go to your repository → Actions tab
   - The "Documentation" workflow will run automatically
   - Wait for it to complete successfully (creates the `gh-pages` branch)

6. **Enable GitHub Pages** (after the first successful workflow run):
   - Go to repository settings → Pages
   - Set source to "Deploy from a branch"
   - Select the `gh-pages` branch and `/ (root)` folder
   - Click Save
   
   **Alternative (recommended)**: Use GitHub Actions source
   - Set source to "GitHub Actions"
   - No need to select a branch - it will deploy automatically

The documentation will automatically build and deploy on every push to the main branch.

## Documentation URL

Once deployed, your documentation will be available at:
https://CarlosP24.github.io/FullShell.jl/dev/

(For tagged releases, it will also be at `https://CarlosP24.github.io/FullShell.jl/stable/`)

## Updating Documentation

### Adding/Modifying Docstrings
1. Edit the docstrings in your Julia source files
2. Follow the format shown in existing docstrings
3. Rebuild documentation to see changes

### Adding Examples
Edit `docs/src/examples.md` to add new examples

### Modifying Pages
Edit the markdown files in `docs/src/` as needed

### Adding New Pages
1. Create a new `.md` file in `docs/src/`
2. Add it to the `pages` array in `docs/make.jl`
3. Rebuild documentation

## Testing Documentation Locally

Before pushing, always test that:
1. Documentation builds without errors
2. All docstrings render correctly
3. Examples are accurate
4. Links work properly

## Tips for Good Documentation

1. **Keep docstrings up-to-date**: Update them when you change function signatures
2. **Provide examples**: Show actual usage with real parameters
3. **Explain physics**: Help users understand what parameters mean physically
4. **Cross-reference**: Link related functions in docstrings
5. **Use LaTeX math**: Wrap equations in backticks for inline or double-dollar signs for display

## For Package Registration

When you're ready to register the package:

1. Make sure documentation is building successfully
2. Verify all docstrings are complete
3. Check that the README links to documentation
4. Follow the Julia package registration process
5. Update the README installation instructions once registered

The current documentation setup is fully compatible with the Julia package registry requirements.

## Troubleshooting

### Build Errors
- Check that all docstrings use proper syntax
- Ensure all `@docs` blocks reference exported functions
- Verify there are no syntax errors in markdown files

### Missing Functions in API
- Make sure the function is exported in the main module
- Check that the docstring is attached to the function definition
- Rebuild with `checkdocs = :exports` to catch issues

### GitHub Actions Failing
- Check the Actions tab in GitHub for error messages
- Verify DOCUMENTER_KEY is set correctly
- Ensure docs/Project.toml has all necessary dependencies

## Resources

- [Documenter.jl Documentation](https://documenter.juliadocs.org/stable/)
- [Julia Documentation Guide](https://docs.julialang.org/en/v1/manual/documentation/)
- [Example packages with good docs](https://github.com/JuliaDocs)
