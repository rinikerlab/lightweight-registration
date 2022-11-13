from setuptools import setup

setup(
    name='lwreg',
    description='lightweight compound registration system',
    version='0.0.1',
    py_modules=['lwreg'],
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'lwreg = lwreg.lwreg:cli',
        ],
    },
    author="Greg Landrum",
    author_email="greg.landrum@gmail.com",

)