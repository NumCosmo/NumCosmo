# NumCosmo Testing Guide

## Test Organization

### Test Directories
- **tests/** - Main test suite directory
  - C tests: `test_*.c` files
  - Python tests: `test_py_*.py` files
  - Test fixtures: `conftest.py`, `fixtures_*.py`

### Test Categories
- Unit tests for individual components
- Integration tests for module interactions
- Regression tests for bug fixes
- Performance benchmarks

## C Testing

### Test Framework
- Uses GLib testing framework (`g_test_*` functions)
- Test registration with `g_test_add_func()`
- Assertions: `g_assert()`, `g_assert_cmpfloat()`, etc.

### Test Structure
```c
#include <numcosmo/numcosmo.h>

static void
test_component_feature (void)
{
  // Setup
  NcmObject *obj = ncm_object_new ();
  
  // Test
  gdouble result = ncm_object_compute (obj);
  
  // Assert
  g_assert_cmpfloat (result, ==, expected_value);
  
  // Cleanup
  ncm_object_free (obj);
}

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  
  g_test_add_func ("/component/feature", test_component_feature);
  
  return g_test_run ();
}
```

### Running C Tests
```bash
# Run all tests
meson test -C builddir

# Run specific test
meson test -C builddir test_name

# Verbose output
meson test -C builddir --verbose

# Run with valgrind
meson test -C builddir --wrap='valgrind --leak-check=full'
```

## Python Testing

### Test Framework
- Uses pytest framework
- Fixtures defined in `conftest.py`
- Parametrized tests with `@pytest.mark.parametrize`

### Test Structure
```python
import pytest
from numcosmo_py import Nc, Ncm

def test_component_feature():
    # Setup
    obj = Ncm.Object.new()
    
    # Test
    result = obj.compute()
    
    # Assert
    assert result == expected_value
    
    # Cleanup (automatic via garbage collection)
```

### Fixtures
- Common fixtures in `conftest.py`
- Specialized fixtures in `fixtures_*.py`
- CCL comparison fixtures in `fixtures_ccl.py`
- Cross-correlation fixtures in `fixtures_xcor.py`

### Running Python Tests
```bash
# Run all Python tests
pytest tests/

# Run specific test file
pytest tests/test_py_hicosmo.py

# Run specific test function
pytest tests/test_py_hicosmo.py::test_hicosmo_de

# Run with coverage
pytest --cov=numcosmo_py tests/

# Run with verbose output
pytest -v tests/
```

## Test Naming Conventions

### C Tests
- File: `test_[module]_[component].c`
- Function: `test_[component]_[feature]()`
- Examples:
  - `test_ncm_spline.c` → `test_spline_interpolation()`
  - `test_nc_distance.c` → `test_distance_comoving()`

### Python Tests
- File: `test_py_[module].py`
- Function: `test_[component]_[feature]()`
- Examples:
  - `test_py_hicosmo.py` → `test_hicosmo_de()`
  - `test_py_powspec.py` → `test_powspec_ml_transfer()`

## Test Coverage

### Generating Coverage Reports
```bash
# Configure with coverage
meson setup -Db_coverage=true builddir

# Run tests
meson test -C builddir

# Generate coverage report
ninja coverage -C builddir

# View HTML report
xdg-open builddir/meson-logs/coveragereport/index.html
```

### Coverage Configuration
- Configuration in `.lcovrc`
- Excludes test files and external libraries
- Target: >80% coverage for core modules

## Continuous Integration

### GitHub Actions
- Workflow: `.github/workflows/build-and-check.yml`
- Runs on: Linux, macOS
- Tests: C and Python test suites
- Coverage: Uploaded to codecov.io

### CI Test Matrix
- Multiple Python versions (3.9, 3.10, 3.11, 3.12, 3.13)
- Different build configurations (debug, release)
- Optional dependencies (with/without FFTW, cfitsio, NLopt)
- MPI builds

## Test Data

### Data Files
- Located in `data/` directory
- Serialized objects: `*.obj` files
- Observational data: BAO, Hubble, SNIa, CMB
- Test-specific data generated during tests

### Using Test Data
```c
// C: Load serialized object
NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
NcmObject *obj = ncm_serialize_from_file (ser, "data/object.obj");
```

```python
# Python: Load serialized object
ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
obj = ser.from_file("data/object.obj")
```

## Debugging Tests

### Debug Build
```bash
# Configure debug build
meson setup --buildtype=debug builddir

# Run test with gdb
gdb --args builddir/tests/test_name
```

### Common Debug Flags
- `G_DEBUG=fatal-warnings` - Make warnings fatal
- `G_SLICE=always-malloc` - Use system malloc
- `G_SLICE=debug-blocks` - Memory debugging

### Valgrind
```bash
# Memory leak check
valgrind --leak-check=full builddir/tests/test_name

# With suppressions
valgrind --suppressions=numcosmo-valgrind.supp builddir/tests/test_name
```

## Best Practices

### Test Design
- Each test should be independent
- Use fixtures for common setup
- Test edge cases and error conditions
- Keep tests fast (< 1 second per test)

### Assertions
- Use appropriate assertion functions
- Include descriptive messages
- Test both success and failure paths

### Cleanup
- Always free allocated resources
- Use `_clear()` macros for GObjects
- Verify no memory leaks with valgrind

### Documentation
- Document test purpose in comments
- Explain non-obvious test logic
- Reference related issues/bugs for regression tests
