from setuptools import setup, find_packages

def calculate_version(inputfile):
    version_list =  [x.split('\'')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

package_version = calculate_version('./pdxSPA/_version.py')

setup(
    name='pdxSPA',
    version=package_version,
    author='Xi Cheng',
    author_email='chengxi0237@sjtu.edu.cn,',
    packages=find_packages(),
    install_requires=[
        'Click',
        'numpy',
        'pandas',
    ],

    entry_points="""
    [console_scripts]
    pdxSPA=pdxSPA.cli:main
    """,

    description='SPA (Shared Peptides Allocation) for quantifying the MS data from the PDX models.',
    license='GPL-3.0',
    classifiers=[
        'Intended Audience :: Proteomics',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    url='https://github.com/Li-Lab-Proteomics/pdxSPA/',

)
