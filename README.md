Phys 581 - The Standard Model of Particle Physics
=================================================
[![Documentation Status](https://readthedocs.org/projects/physics-581-the-standard-model/badge/?version=latest)](https://physics-581-the-standard-model.readthedocs.io/en/latest/?badge=latest)
[![gitlab pipeline status](https://gitlab.com/wsu-courses/physics-581-the-standard-model/badges/main/pipeline.svg)](https://gitlab.com/wsu-courses/physics-581-the-standard-model/-/commits/main)
[![gitlab coverage report](https://gitlab.com/wsu-courses/physics-581-the-standard-model/badges/main/coverage.svg)](https://gitlab.com/wsu-courses/physics-581-the-standard-model/-/commits/main)

This is the main project for the [WSU Physics][] course
**Phys 581: The Standard Model of Particle Physics** first offered in [Fall 2023](https://schedules.wsu.edu/List/Pullman/20233/Phys/581/01).

Physics has a successful track record of providing effective solutions to complex
problems outside its specific domain. This course will focus on using efficient
numerical techniques inspired by physics to solve challenging problems in a wide variety
of applications.  Techniques will be chosen from physics applications, but also applied
to problems outside of the physics domain including economics, biology, sociology, etc.
Students will be introduced to powerful numerical toolkits based on the
[SciPy](https://www.scipy.org/) and [NumFocus](https://numfocus.org) ecosystem. Using
the [CoCalc](https://cocalc.com/) platform will enable rapid development and prototyping
with an explicit path to stable, tested, and performant codes capable of supporting
research, or industry applications.

## Docs

To build the documents interactively:

```bash
make doc-server
```

This will run [`sphinx-autobuild`](https://github.com/executablebooks/sphinx-autobuild)
which will launch a webserver on http://127.0.0.1:8000 and rebuild the docs whenever you
save a change.

Here is the play-by-play for setting up the documentation.

```bash
cd Docs
sphinx-quickstart
wget https://brand.wsu.edu/wp-content/themes/brand/images/pages/logos/wsu-signature-vertical.svg -O _static/wsu-logo.svg 
cp -r ../envs/default/lib/python3.9/site-packages/sphinx_book_theme/_templates/* _templates
```

I then edited the `conf.py`

```bash
hg add local.bib _static/ _templates/
```

## CoCalc Setup


* [Purchase a license](https://cocalc.com/settings/licenses) with 2 projects to allow
  the course and [WSU Courses CoCalc project][] and [Shared CoCalc Project][] to run.  This
  approach requires the students to pay $14 for access four the term (4 months).  They
  can optionally use any license they already have instead.
   
  Optionally, one might opt to purchase a license for $n+2$ projects where $n$ is the
  number of students, if there is central funding available.  See [Course Upgrading
  Students](https://doc.cocalc.com/teaching-upgrade-course.html#course-upgrading-students)
  for more details.
  
* Next, [create a course](https://doc.cocalc.com/teaching-create-course.html).  I do
  this in my [WSU Courses CoCalc project][].



* Create a [Shared CoCalc Project][] and activate the license for this project so that it
  can run.  I then add the SSH key to may `.ssh/config` files so I can quickly login.

* Clone the repos into the shared project and initialize the project.  Optional, but
  highly recommend -- use my [`mmf-setup`][] project to provide some useful features

  ```bash
  ssh smc<project name>       # My alias in .ssh/config
  python3 -m pip install mmf_setup
  mmf_setup cocalc
  ```
  
  This provides some instructions on how to use the CoCalc configuration.  The most
  important is to forward your user agent and set your `hg` and `git` usernames:
  
  ```bash
  ~$ mmf_setup cocalc
  ...
  If you use version control, then to get the most of the configuration,
  please make sure that you set the following variables on your personal
  computer, and forward them when you ssh to the project:

      # ~/.bashrc or similar
      LC_HG_USERNAME=Your Full Name <your.email.address+hg@gmail.com>
      LC_GIT_USEREMAIL=your.email.address+git@gmail.com
      LC_GIT_USERNAME=Your Full Name

  To forward these, your SSH config file (~/.ssh/config) might look like:

      # ~/.ssh/config
      Host cc_phys-581-the-standard-model-of-particle-physics
        User 4578ccc855cf413d9e4c82c69e30121e
    
      Host cc_phys-581-the-standard-model-of-particle-physics
        HostName ssh.cocalc.com
        ForwardAgent yes
        SendEnv LC_HG_USERNAME
        SendEnv LC_GIT_USERNAME
        SendEnv LC_GIT_USEREMAIL
        SetEnv LC_EDITOR=vi
  ```
  
  Logout and log back in so we have the forwarded credentials, and now clone the repos.
  
  ```bash
  git clone https://gitlab.com/wsu-courses/physics-581-the-standard-model.git phys-581-the-standard-model-of-particle-physics
  cd phys-581-the-standard-model-of-particle-physics
  make
  ```
  
  The last step runs `git clone git@gitlab.com:wsu-courses/physics-581-the-standard-model_resources.git _ext/Resources` which puts the resources folder in `_ext/Resources`.

* Create an environment:

  ```bash
  ssh cc_phys-581-the-standard-model-of-particle-physics
  cd phys-581-the-standard-model-of-particle-physics
  anaconda2021
  anaconda-project prepare
  conda activate envs/phys-581
  python -m ipykernel install --user --name "phys-581" --display-name "Python 3 (phys-581)"
  ```

  This will create a Conda environment as specified in `anaconda-project.yml` in `envs/phys-581`.


# Funding

<a href="https://www.nsf.gov"><img width="10%"
src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png" />
</a>
<br>

Some of the material presented here is based upon work supported by the National Science
Foundation under [Grant Number 2012190](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2012190). Any opinions, findings, and conclusions or
recommendations expressed in this material are those of the author(s) and do not
necessarily reflect the views of the National Science Foundation.


<!-- Links -->
[Anaconda Project]: <https://github.com/Anaconda-Platform/anaconda-project> "Anaconda Project"
[CoCalc]: <https://cocalc.com> "CoCalc: Collaborative Calculation and Data Science"
[Conda]: <https://docs.conda.io/en/latest/> "Conda: Package, dependency and environment management for any language—Python, R, Ruby, Lua, Scala, Java, JavaScript, C/ C++, FORTRAN, and more."
[GitHub CI]: <https://docs.github.com/en/actions/guides/about-continuous-integration> "GitHub CI"
[GitHub]: <https://github.com> "GitHub"
[GitLab]: <https://gitlab.com> "GitLab"
[Git]: <https://git-scm.com> "Git"
[Heptapod]: <https://heptapod.net> "Heptapod: is a community driven effort to bring Mercurial SCM support to GitLab"
[Jupyter]: <https://jupyter.org> "Jupyter"
[Jupytext]: <https://jupytext.readthedocs.io> "Jupyter Notebooks as Markdown Documents, Julia, Python or R Scripts"
[Mercurial]: <https://www.mercurial-scm.org> "Mercurial"
[Miniconda]: <https://docs.conda.io/en/latest/miniconda.html> "Miniconda is a free minimal installer for conda."
[MyST]: <https://myst-parser.readthedocs.io/en/latest/> "MyST - Markedly Structured Text"
[Read the Docs]: <https://readthedocs.org> "Read the Docs homepage"
[WSU Physics]: <https://physics.wsu.edu> "WSU Department of Physics and Astronomy"
[WSU Mathematics]: <https://www.math.wsu.edu/> "WSU Department of Mathematics and Statistics"
[`anaconda-project`]: <https://anaconda-project.readthedocs.io> "Anaconda Project: Tool for encapsulating, running, and reproducing data science projects."
[`anybadge`]: <https://github.com/jongracecox/anybadge> "Python project for generating badges for your projects"
[`conda-forge`]: <https://conda-forge.org/> "A community-led collection of recipes, build infrastructure and distributions for the conda package manager."
[`genbadge`]: <https://smarie.github.io/python-genbadge/> "Generate badges for tools that do not provide one."
[`mmf-setup`]: <https://pypi.org/project/mmf-setup/> "PyPI mmf-setup page"
[`pytest`]: <https://docs.pytest.org> "pytest: helps you write better programs"
[hg-git]: <https://hg-git.github.io> "The Hg-Git mercurial plugin"
[GitLab test coverage visualization]: <https://docs.gitlab.com/ee/user/project/merge_requests/test_coverage_visualization.html>



[Official Course Repository]: <https://gitlab.com/wsu-courses/physics-581-the-standard-model> "Official Course Repository hosted on GitLab"
[Shared CoCalc Project]: <https://cocalc.com/4578ccc8-55cf-413d-9e4c-82c69e30121e/> "Shared CoCalc Project"
[WSU Courses CoCalc project]: <https://cocalc.com/projects/c31d20a3-b0af-4bf7-a951-aa93a64395f6>


[GitHub Mirror]: <https://github.com/WSU-Physics-Courses/physics-581-the-standard-model> "GitHub mirror"
[GitLab public repo]: <https://gitlab.com/wsu-courses/physics-581-the-standard-model> "GitLab public repository."
[Gitlab private resources repo]: <https://gitlab.com/wsu-courses/physics-581-the-standard-model_resources> "Private resources repository."
[file an issue]: <https://gitlab.com/wsu-courses/physics-581-the-standard-model/-/issues> "Issues on the GitLab project."

<!-- End Links -->
