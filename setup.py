from setuptools import setup, find_packages

setup(
    name="EBM1DModel",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'pandas==1.4.2',
        'numpy==1.22.3',
        'matplotlib==3.5.2',
        'mpmath==1.2.1'
    ],
    entry_points={
        'console_scripts': [
            'run_ebm=EBM_1D:main',
        ],
    },
    package_data={
        '': ['*.dat'],
    },
    include_package_data=True,
    long_description="""
    A 1-dimensional energy balance model adapted from the work found at:
    https://github.com/missrg11/1D-Energy-Balance-Model
    
    This model simulates the surface temperature of Earth across different latitude zones,
    accounting for variations in solar energy, albedo effects, and other climate-related factors.
    """
)
