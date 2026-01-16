# Contributing to FullShell.jl

Thank you for considering contributing to FullShell.jl! This document provides guidelines for contributing to the project.

## Development Setup

1. Clone the repository:
```bash
git clone https://github.com/CarlosP24/FullShell.jl.git
cd FullShell.jl
```

2. Install the package in development mode:
```julia
using Pkg
Pkg.develop(path=".")
```

3. Run tests:
```julia
Pkg.test("FullShell")
```

## Building Documentation Locally

To build and preview the documentation locally:

1. Navigate to the docs directory:
```bash
cd docs
```

2. Instantiate the docs environment:
```julia
using Pkg
Pkg.activate(".")
Pkg.develop(PackageSpec(path=".."))
Pkg.instantiate()
```

3. Build the documentation:
```julia
include("make.jl")
```

4. The built documentation will be in `docs/build/`. Open `docs/build/index.html` in your browser to preview.

Alternatively, you can use a local server:
```bash
# From the docs/build directory
python -m http.server 8000
# Then open http://localhost:8000 in your browser
```

## Adding Documentation

### Docstrings

All exported functions should have docstrings following the Julia documentation standard:

```julia
"""
    function_name(arg1, arg2; kwarg=default)

Brief description of what the function does.

# Arguments
- `arg1::Type`: Description of arg1
- `arg2::Type`: Description of arg2
- `kwarg::Type`: Description of keyword argument (default: `default`)

# Returns
- Description of return value(s)

# Examples
\```julia
result = function_name(1, 2)
\```
"""
function function_name(arg1, arg2; kwarg=default)
    # implementation
end
```

### Documentation Pages

Documentation pages are located in `docs/src/` and written in Markdown with Documenter.jl extensions:

- `index.md`: Main landing page
- `examples.md`: Usage examples
- `api.md`: Complete API reference

To add examples, edit `docs/src/examples.md`. To add new API documentation, ensure your functions have proper docstrings and are exported from the main module.

## Code Style

- Follow [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- Use meaningful variable names
- Add comments for complex algorithms
- Keep functions focused and concise

## Testing

- Add tests for new functionality in the `test/` directory
- Ensure all tests pass before submitting a pull request
- Aim for good test coverage of new code

## Pull Request Process

1. Fork the repository and create a new branch for your feature
2. Make your changes, including tests and documentation
3. Ensure all tests pass
4. Update CHANGELOG.md (if applicable)
5. Submit a pull request with a clear description of the changes

## Reporting Issues

When reporting issues, please include:
- Julia version
- FullShell.jl version
- Minimal reproducible example
- Expected vs actual behavior
- Any error messages or stack traces

## Questions

If you have questions, feel free to:
- Open an issue for discussion
- Contact the maintainers

Thank you for contributing!
