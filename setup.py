from setuptools import setup

setup(name='densityplot',
      version='1.0.0',
      description='tools to plot large numbers of scatter points',
      author='Dr. Coleman Krawczyk',
      author_email='coleman.krawczyk@gmail.com',
      packages=['densityplot'],
      install_requires=['matplotlib','scipy','numpy'],
      zip_safe=False
)
