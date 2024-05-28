from setuptools import setup, find_packages

setup(
    name='scRNA_clusterings',
    version='0.1',
    description='A package for single-cell RNA sequencing analysis with various clustering methods.',
    author='Shiyuan Wen',
    author_email='swen@ucsd.edu',
    url='https://github.com/ShiyuanWen/CSE185_FinalProject/tree/main',
    packages=find_packages(),
    install_requires=[
        'scanpy',
        'anndata',
        'harmonypy',
        'matplotlib',
        'seaborn',
        'leidenalg',
        'sc3s',
        'scikit-learn'
    ],
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
