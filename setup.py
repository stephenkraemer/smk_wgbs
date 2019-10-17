from setuptools import setup

setup(
        name='scwgbs_alignment',
        version='0.1',
        author='Stephen Kraemer',
        author_email='stephenkraemer@gmail.com',
        license='MIT',
        packages=['smk_wgbs'],
        package_data={
            '': ['*.R', '*.snakefile', '*.yml', '*.yaml', '*.sh'],
        },
        install_requires=[],
        python_requires='>=3.6',
)
