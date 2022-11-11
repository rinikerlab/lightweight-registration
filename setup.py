from setuptools import setup

setup(
    name='lwreg',
    version='0.0.1',
    py_modules=['lwreg'],
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'lwreg = lwreg:cli',
        ],
    },
)