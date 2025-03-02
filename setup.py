from setuptools import setup, find_packages

setup(
    name="ribxz-chimerax",  # Use hyphens instead of underscores
    version="0.1",
    package_dir={"": "src"},  # Tell setuptools packages are under src/
    packages=find_packages(where="src"),
    package_data={
        "ribxz_chimerax": ["bundle_info.xml"],
    },
    install_requires=[
        "aiohttp",
    ],
    entry_points={
        "chimerax.bundle": ["ribxz_chimerax = ribxz_chimerax:bundle_api"],
    },
)