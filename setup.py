from setuptools import setup, find_packages

VERSION = '0.0.7' 
DESCRIPTION = 'Numerical simulation of some PDEs using finite difference methods'
LONG_DESCRIPTION = 'pic16b'

# Setting up
setup(
       
        name="PdeSimulation", 
        version=VERSION,
        author="Yun",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'first package'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)
