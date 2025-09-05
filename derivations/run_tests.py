#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Physics Verification Test Runner
================================

Automated test runner for the physics verification framework. Discovers and runs
all standardized test modules, providing comprehensive validation of the mathematical
correctness of the vortex field theory framework.

This runner assumes all test modules follow the standard structure defined in
TEST_STANDARD.md, with each module having a test_[module_name]() function that
returns the success rate.

Usage:
    python run_tests.py          # Run all tests with normal output
    python run_tests.py --quiet  # Run tests in quiet mode (summary only)
    python run_tests.py --help   # Show help information
"""

import os
import sys
import importlib.util
import argparse
import traceback
import time
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any

# Try to import helper for global test counting
try:
    from helper import reset_global_test_counts, get_global_test_counts
    _HELPER_AVAILABLE = True
except ImportError:
    _HELPER_AVAILABLE = False


class TestResult:
    """Container for individual test module results."""

    def __init__(self, module_name: str, file_path: str):
        self.module_name = module_name
        self.file_path = file_path
        self.success_rate = 0.0
        self.execution_time = 0.0
        self.error = None
        self.completed = False
        self.tests_passed = 0
        self.tests_total = 0

    def set_success(self, success_rate: float, execution_time: float, tests_passed: int = 0, tests_total: int = 0):
        """Mark test as successfully completed."""
        self.success_rate = success_rate
        self.execution_time = execution_time
        self.completed = True
        self.tests_passed = tests_passed
        self.tests_total = tests_total

    def set_error(self, error_message: str):
        """Mark test as failed with error."""
        self.error = error_message
        self.completed = False

    def is_success(self) -> bool:
        """Check if test completed successfully with 100% pass rate."""
        return self.completed and self.success_rate == 100.0

    def is_completed(self) -> bool:
        """Check if test completed (regardless of pass rate)."""
        return self.completed

    def __repr__(self):
        if self.error:
            return f"TestResult({self.module_name}: ERROR)"
        elif self.completed:
            return f"TestResult({self.module_name}: {self.success_rate:.1f}%)"
        else:
            return f"TestResult({self.module_name}: NOT RUN)"


class PhysicsTestRunner:
    """
    Main test runner for standardized physics verification tests.

    Discovers test modules, runs their test_[module_name]() functions,
    and aggregates results into a comprehensive report.
    """

    def __init__(self, quiet: bool = False, verbose: bool = False):
        self.quiet = quiet
        self.verbose = verbose
        self.results: List[TestResult] = []

        # Directories to scan for tests (excluding emergent_particle_masses as requested)
        self.test_directories = [
            "mathematical_framework",
            "projected_em",
            "emergent_particle_masses",
            "gravity",
            "quantum"
        ]

        # Files to exclude from testing
        self.excluded_files = {
            "helper.py",
            "helper.test.py",
            "__init__.py",
            "run_tests.py",
            "TEST_STANDARD.md"
        }

    def discover_test_modules(self) -> List[str]:
        """
        Discover all Python test modules in the specified directories.

        Returns:
            List of module file paths relative to current directory
        """
        test_files = []
        # Use the directory where this script is located, not the current working directory
        script_dir = Path(__file__).parent
        base_dir = script_dir

        for test_dir in self.test_directories:
            dir_path = base_dir / test_dir
            if dir_path.exists() and dir_path.is_dir():
                for py_file in dir_path.glob("*.py"):
                    if py_file.name not in self.excluded_files:
                        test_files.append(str(py_file))
                        if self.verbose:
                            print(f"Discovered test module: {py_file}")

        return sorted(test_files)

    def import_test_module(self, module_path: str) -> Optional[Any]:
        """
        Dynamically import a test module from its file path.

        Args:
            module_path: Path to the Python module file

        Returns:
            Imported module or None if import failed
        """
        try:
            # Generate module name from path
            module_name = Path(module_path).stem
            spec = importlib.util.spec_from_file_location(module_name, module_path)

            if spec is None or spec.loader is None:
                raise ImportError(f"Could not create module spec for {module_path}")

            module = importlib.util.module_from_spec(spec)

            # Ensure script directory is in Python path for helper imports
            script_dir = str(Path(__file__).parent.resolve())
            if script_dir not in sys.path:
                sys.path.insert(0, script_dir)

            spec.loader.exec_module(module)
            return module

        except Exception as e:
            if self.verbose:
                print(f"Import error for {module_path}: {e}")
                print(traceback.format_exc())
            return None

    def run_single_test(self, module_path: str) -> TestResult:
        """
        Run a single test module following the standard structure.

        Args:
            module_path: Path to the test module

        Returns:
            TestResult with test outcomes
        """
        module_name = Path(module_path).stem
        result = TestResult(module_name, module_path)

        if not self.quiet:
            print(f"\n{'='*60}")
            print(f"Running: {module_name}")
            print('='*60)
        elif self.verbose:
            print(f"Running {module_name}...")

        try:
            # Import the module
            start_time = time.time()
            module = self.import_test_module(module_path)

            if module is None:
                result.set_error("Failed to import module")
                return result

            # Look for the standard test function
            test_function_name = f"test_{module_name}"

            if not hasattr(module, test_function_name):
                result.set_error(f"Module missing required function: {test_function_name}()")
                return result

            test_function = getattr(module, test_function_name)

            if not callable(test_function):
                result.set_error(f"{test_function_name} is not callable")
                return result

            # Set quiet flag in sys.argv for helper.py to detect
            quiet_added = False
            if self.quiet and '--quiet' not in sys.argv:
                sys.argv.append('--quiet')
                quiet_added = True

            try:
                # Reset global test counters before running the test
                if _HELPER_AVAILABLE:
                    reset_global_test_counts()

                # Run the test function
                success_rate = test_function()
                execution_time = time.time() - start_time

                # Get global test counts after running the test
                tests_passed, tests_total = 0, 0
                if _HELPER_AVAILABLE:
                    tests_passed, tests_total = get_global_test_counts()

                # Validate return value
                if not isinstance(success_rate, (int, float)):
                    result.set_error(f"Test function returned {type(success_rate)}, expected numeric success rate")
                    return result

                if not (0 <= success_rate <= 100):
                    result.set_error(f"Success rate {success_rate} outside valid range [0, 100]")
                    return result

                result.set_success(float(success_rate), execution_time, tests_passed, tests_total)

            finally:
                # Clean up sys.argv
                if quiet_added:
                    sys.argv.remove('--quiet')

        except Exception as e:
            execution_time = time.time() - start_time
            error_msg = f"Exception during test execution: {e}"
            if self.verbose:
                error_msg += f"\n{traceback.format_exc()}"
            result.set_error(error_msg)

        return result

    def run_all_tests(self) -> Dict[str, Any]:
        """
        Run all discovered tests and return comprehensive results.

        Returns:
            Dictionary with detailed test results and statistics
        """
        print("Physics Verification Framework - Standardized Test Runner")
        print("="*60)

        if self.quiet:
            print("Running in quiet mode - showing summary only\n")

        # Discover test modules
        test_modules = self.discover_test_modules()

        if not test_modules:
            print("No test modules found!")
            return {
                "total_modules": 0,
                "completed_modules": 0,
                "successful_modules": 0,
                "overall_success_rate": 0.0,
                "results": []
            }

        print(f"Discovered {len(test_modules)} test modules")
        if not self.quiet:
            for module in test_modules:
                print(f"  - {module}")
        print()

        # Run each test module
        for module_path in test_modules:
            result = self.run_single_test(module_path)
            self.results.append(result)

            # Show immediate feedback in quiet mode
            if self.quiet:
                status = "✅" if result.is_success() else "❌" if result.is_completed() else "💥"
                rate = f"{result.success_rate:.1f}%" if result.is_completed() else "ERROR"
                test_info = f" ({result.tests_passed}/{result.tests_total})" if result.tests_total > 0 else ""
                print(f"{status} {result.module_name:<40} {rate}{test_info}")

        # Generate and return comprehensive summary
        return self.generate_summary()

    def generate_summary(self) -> Dict[str, Any]:
        """Generate comprehensive summary of all test results."""
        total_modules = len(self.results)
        completed_modules = sum(1 for r in self.results if r.is_completed())
        successful_modules = sum(1 for r in self.results if r.is_success())
        failed_modules = [r for r in self.results if not r.is_completed()]

        # Calculate weighted overall success rate (only from completed modules)
        if completed_modules > 0:
            overall_success_rate = sum(r.success_rate for r in self.results if r.is_completed()) / completed_modules
        else:
            overall_success_rate = 0.0

        total_execution_time = sum(r.execution_time for r in self.results)

        # Calculate total individual test counts
        total_individual_tests = sum(r.tests_total for r in self.results if r.is_completed())
        total_individual_passed = sum(r.tests_passed for r in self.results if r.is_completed())

        summary = {
            "total_modules": total_modules,
            "completed_modules": completed_modules,
            "successful_modules": successful_modules,
            "failed_modules": len(failed_modules),
            "overall_success_rate": overall_success_rate,
            "module_success_rate": (successful_modules / total_modules * 100) if total_modules > 0 else 0.0,
            "total_execution_time": total_execution_time,
            "total_individual_tests": total_individual_tests,
            "total_individual_passed": total_individual_passed,
            "individual_success_rate": (total_individual_passed / total_individual_tests * 100) if total_individual_tests > 0 else 0.0,
            "results": self.results,
            "failed_results": failed_modules
        }

        self.print_summary(summary)
        return summary

    def print_summary(self, summary: Dict[str, Any]):
        """Print formatted comprehensive summary of all test results."""
        print("\n" + "="*80)
        print("COMPREHENSIVE TEST SUMMARY")
        print("="*80)

        # Overall statistics
        print(f"Total test modules: {summary['total_modules']}")
        print(f"Modules completed: {summary['completed_modules']}")
        print(f"Modules with 100% pass rate: {summary['successful_modules']}")
        print(f"Modules with errors: {summary['failed_modules']}")
        print(f"Overall success rate: {summary['overall_success_rate']:.1f}%")
        print(f"Module success rate: {summary['module_success_rate']:.1f}%")
        if summary['total_individual_tests'] > 0:
            print(f"Individual tests run: {summary['total_individual_tests']}")
            print(f"Individual tests passed: {summary['total_individual_passed']}")
            print(f"Individual test success rate: {summary['individual_success_rate']:.1f}%")
        print(f"Total execution time: {summary['total_execution_time']:.2f}s")

        # Per-module results (skip in quiet mode)
        if not self.quiet:
            print(f"\nDETAILED RESULTS:")
            print("-" * 70)
            print(f"{'STATUS':<8} {'MODULE':<35} {'SUCCESS RATE':<12} {'TIME':<8}")
            print("-" * 70)

            for result in self.results:
                if result.error:
                    status = "❌ ERROR"
                    rate = "N/A"
                    time_str = "N/A"
                elif result.is_success():
                    status = "✅ PASS"
                    rate = f"{result.success_rate:.1f}%"
                    time_str = f"{result.execution_time:.2f}s"
                elif result.is_completed():
                    status = "⚠️ PARTIAL"
                    rate = f"{result.success_rate:.1f}%"
                    time_str = f"{result.execution_time:.2f}s"
                else:
                    status = "💥 FAILED"
                    rate = "N/A"
                    time_str = "N/A"

                print(f"{status:<8} {result.module_name:<35} {rate:<12} {time_str:<8}")

        # Show error details if any
        if summary['failed_results']:
            print(f"\nERROR DETAILS:")
            print("-" * 50)
            for result in summary['failed_results']:
                print(f"❌ {result.module_name}:")
                print(f"   {result.error}")
                if result.file_path:
                    print(f"   File: {result.file_path}")
                print()

        # Overall assessment (skip in quiet mode)
        if not self.quiet:
            print(f"OVERALL ASSESSMENT:")
            print("-" * 40)

            if summary['failed_modules'] > 0:
                print("❌ SOME MODULES FAILED TO RUN")
                print(f"   {summary['failed_modules']} modules had import or execution errors")
            elif summary['overall_success_rate'] == 100.0:
                print("🎉 ALL VERIFICATIONS PASSED!")
                print("   Mathematical framework is fully verified.")
            elif summary['overall_success_rate'] >= 95.0:
                print("✅ VERIFICATION SUBSTANTIALLY COMPLETE")
                print(f"   Average success rate: {summary['overall_success_rate']:.1f}%")
            elif summary['overall_success_rate'] >= 85.0:
                print("⚠️  VERIFICATION MOSTLY SUCCESSFUL")
                print(f"   Average success rate: {summary['overall_success_rate']:.1f}%")
            else:
                print("❌ VERIFICATION NEEDS ATTENTION")
                print(f"   Average success rate: {summary['overall_success_rate']:.1f}%")

            print("\nThis test runner validates the mathematical correctness of")
            print("the vortex field theory framework documented in ../doc/")
            print("\nAll tests follow the standardized structure defined in TEST_STANDARD.md")


def main():
    """Main entry point for the standardized test runner."""
    parser = argparse.ArgumentParser(
        description="Run all standardized physics verification tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This test runner assumes all test modules follow the standard structure
defined in TEST_STANDARD.md, with each module having a main function
test_[module_name]() that returns a success rate (0-100).

Examples:
    python run_tests.py              # Normal output with detailed results
    python run_tests.py --quiet      # Quiet mode - summary only
    python run_tests.py --verbose    # Verbose mode with debug information

The runner scans mathematical_framework/ and projected_em/ directories
for Python files and runs their standardized test functions.
        """
    )

    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Run in quiet mode (minimal output, summary only)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Run in verbose mode with detailed debug information')

    args = parser.parse_args()

    if args.quiet and args.verbose:
        print("Error: --quiet and --verbose flags are mutually exclusive")
        sys.exit(1)

    # Create and run test runner
    runner = PhysicsTestRunner(quiet=args.quiet, verbose=args.verbose)
    summary = runner.run_all_tests()

    # Set exit code based on results
    # Exit with failure if any modules failed to run or average success rate < 100%
    if summary['failed_modules'] > 0 or summary['overall_success_rate'] < 100.0:
        sys.exit(1)  # Failure exit code
    else:
        sys.exit(0)  # Success exit code


if __name__ == "__main__":
    main()
