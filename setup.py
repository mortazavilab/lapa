from setuptools import setup, find_packages


requirements = [
    'setuptools',
    'tqdm',
    'click',
    'pyranges',
    'pandas',
    'kipoiseq',
    'pybigwig',
    'scipy'
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest']

setup(
    author="M. Hasan Çelik",
    author_email='muhammedhasancelik@gmail.com',
    classifiers=[
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
    ],
    description="Fix UTR of Gencode",
    install_requires=requirements,
    license="MIT license",
    entry_points='''
        [console_scripts]
        longread_postprocessing=longread_postprocessing.main:cli
    ''',
    keywords=['genomics', 'long read RNA-seq', 'APA'],
    name='longread_postprocessing',
    packages=find_packages(include=['longread_postprocessing']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/MuhammedHasan/longread_postprocessing',
    version='1.0.0',
)
