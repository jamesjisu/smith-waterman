import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='smith-waterman',
    version='0.0.1',
    author='James Han',
    author_email='jamesjisu@gmail.com',
    description='Smith-Waterman Implementation',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/jamesjisu/smith-waterman',
    license='MIT',
    packages=['smith-waterman'],
    install_requires=['numpy', 'pandas'],
)