from setuptools import setup, find_packages

setup(
    name="quantSim",
    version="0.1.0",
    description="Classical Simulator of Quantum Circuits",
    url="https://https://github.com/sw478/qSim",
    author="Silas Wong",
    author_email="silas.wong478@gmail.com",
    license="MIT",
    classifiers=[
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: OS Independent"
    ],
    packages=["quantSim"],
    include_package_data=True,
    install_requires=["numpy", "dataclasses", "pytest"],
)