from setuptools import setup

setup(
    name='pygsw',
    version='0.1',
    description='Python wrapper for GSW-Fortran',
    author='Jorn Bruggeman',
    author_email='jorn@bolding-bruggeman.com',
    license='GPL',
    packages=['pygsw'],
    package_data={'pygsw': ['*.so', '*.dll', '*.dylib', '*.pyd']},
    zip_safe=False
)


