# Modeled somewhat after
# https://github.com/simoninireland/introduction-to-epidemics/blob/master/Makefile
SHELL = /bin/bash
_SHELL = $(notdir $(SHELL))
DOCS ?= Docs

LOCAL_TOOLS = true

COOKIECUTTER_URL ?= git+https://gitlab.com/forbes-group/cookiecutters.git
COOKIECUTTER ?= cookiecutter $(COOKIECUTTER_URL) --directory project
COOKIECUTTER_YAML ?= .cookiecutter.yaml

GET_RESOURCES = git clone git@gitlab.com:wsu-courses/physics-581-the-standard-model_resources.git _ext/Resources

ENV ?= phys-581
ENV_PATH ?= $(abspath envs/$(ENV))
# Preparation for use on CoCalc.

ifdef ANACONDA2022
	ANACONDA_CURRENT = $(ANACONDA2022)
else ifdef ANACONDA2021
	ANACONDA_CURRENT = $(ANACONDA2021)
else ifdef ANACONDA2020
	ANACONDA_CURRENT = $(ANACONDA2020)
endif

ifdef ANACONDA_CURRENT
	ON_COCALC = true
  COCALC_OPTION ?= micromamba
  COCALC_OPTION ?= anaconda
endif

ifdef ON_COCALC
  ifeq ($(COCALC_OPTION), anaconda)
    # Old approach using anaconda-project in the ANACONDA_CURRENT environment.
    # Due to the $(ANACONDA_CURRENT)/.condarc issue, we must use mamba in this case
    # https://github.com/Anaconda-Platform/anaconda-project/issues/334#issuecomment-911918761
    MAMBA_EXE ?= $(ANACONDA_CURRENT)/bin/mamba
    CREATE_ENV ?= $(CONDA_PRE) $(CONDA_EXE) env create -f $< -p $@
    UPDATE_ENV ?= $(CONDA_PRE) $(CONDA_EXE) env update -f $< -p $@
    ACTIVATE_ENV ?= source $(ANACONDA_CURRENT)/bin/activate && $(CONDA_EXEC) activate -p $(ENV_PATH)
    ANACONDA_PROJECT_PRE ?= CONDA_EXE=$(CONDA_EXE) CONDA_ROOT=$(ANACONDA_CURRENT)
  else
    # New approach - use our own installation of micromamba
    MICROMAMBA_EXE = $(BIN)/micromamba
    MAMBA_EXE = $(MICROMAMBA_EXE)
    CONDA_PRE ?= CONDA_EXE=$(MAMBA_EXE) CONDA_ROOT=~/.local/

    #export MAMBA_ROOT_PREFIX = ~/micromamba
    #MAMBA_PRE = MAMBA_ROOT_PREFIX=$(MAMBA_ROOT_PREFIX)
    _MICROMAMBA = eval "$$($(MAMBA_EXE) shell hook --shell=$(_SHELL) 2> /dev/null)" && \
                  $(CONDA_PRE) $(MAMBA_EXE)
    CREATE_ENV ?= $(_MICROMAMBA) create -f $< -p $@
    UPDATE_ENV ?= $(_MICROMAMBA) update -f $< -p $@
    ACTIVATE_ENV ?= $(_MICROMAMBA) activate $(ENV_PATH)
    INIT_TOOLS += $(MAMBA_EXE)
  endif
  CONDA_EXE ?= $(MAMBA_EXE)
endif

######################################################################
# Generic commands: these use values defined above/
# RUN: run a command like sphinx-build in the appropriate environment

MICROMAMBA_EXE ?= micromamba
MAMBA_EXE = $(MICROMAMBA_EXE)
_MICROMAMBA = eval "$$($(MAMBA_EXE) shell hook --shell=$(_SHELL) 2> /dev/null)" && \
              $(CONDA_PRE) $(MAMBA_EXE)
CREATE_ENV ?= $(_MICROMAMBA) create -y -f $< -p $@

# The following allows micromamba to be run, even if the shell is not initialized.
_RUN = $(_MICROMAMBA) run -p $(ENV_PATH)

CONDA_EXE ?= conda
CONDA ?= eval "$$($(CONDA_PRE) $(CONDA_EXE) shell.$(_SHELL) hook)" && \
         $(CONDA_PRE) $(CONDA_EXE)
ACTIVATE_ENV ?= $(CONDA) activate -p $(ENV_PATH)

INIT_DEPS = ~/.local/bin/mmf_setup

INIT_DEPS += _ext/Resources

ifeq ($(COCALC_OPTION), micromamba)
INIT_DEPS += ~/.local/bin/micromamba
endif
INIT_DEPS += environment.yaml
INIT_DEPS += pyproject.toml
INIT_DEPS += $(ENV_PATH)

ifdef CONDA_SUBDIR
  CONDA_PRE += CONDA_SUBDIR=$(CONDA_SUBDIR)
endif

# ------- Top-level targets  -------
# Default prints a help message
help:
	@make usage

usage:
	@echo "$$HELP_MESSAGE"

info:
	$(MICROMAMBA) info -p $(ENV_PATH)

.PHONY: help usage info

SHELL_INIT_FILE ?= .init-file.$(_SHELL)
shell: $(ENV_PATH)
	$(_RUN) $(SHELL) --init-file $(SHELL_INIT_FILE)

.PHONY: shell

init: $(INIT_TOOLS) $(INIT_DEPS)
ifdef COCALC_ANACONDA
	if ! grep -Fq '$(ACTIVATE_ENV)' ~/.bash_aliases; then \
	  echo '$(ACTIVATE_ENV)' >> ~/.bash_aliases; \
	fi
	make sync
endif
	$(_RUN) python3 -m ipykernel install --user --name "phys-581" \
                                      --display-name "Python 3 (phys-581)"

$(ENV_PATH): environment.yaml pyproject.toml
	$(CREATE_ENV)
	# Would be nice if we can do this in the environment.yaml file.
	$(_RUN) python3 -m pip install --upgrade .[test]

# Jupytext
sync:
	find . -name ".ipynb_checkpoints" -prune -o \
	       -name "_ext" -prune -o \
	       -name "envs" -prune -o \
	       -name "*.ipynb" -o -name "*.md" \
	       -exec $(SHELL) -c '\
               $(JUPYTEXT) --sync "$$1" 2> >(grep -v "is not a paired notebook" 1>&2)' \
               $(SHELL) {} +
# See https://stackoverflow.com/a/15936384/1088938 for details

cleanspace:
	-find . -name "__pycache__" -exec $(RM) -r {} +
	-$(RM) -r _htmlcov .coverage .pytest_cache build
	-$(_MICROMAMBA) clean --all --yes

cleandocs:
	-$(RM) -r $(DOCS)/_build

clean: cleanspace cleandocs

realclean: clean
	$(RM) -r envs

.PHONY: init lock sync clean cleanspace realclean tools

test: $(ENV_PATH)
	$(_RUN) pytest

.PHONY: test

######################################################################
# Tools
#
# This section contains targes for installing tools like micromamba, condax, and
# anaconda-project if needed.  These can either be installed "globally" in
# BIN=~/.local/bin or "locally" if BIN=.local/bin (or anything other than ~/.local/bin).
# The behaviour is as follows:
#
# * On CoCalc, the default is to use BIN=~/.local/bin so that tools can be shared.
#   These WILL be installed with `make init` etc.
# * On other platforms, the default will be to use BIN=.local/bin and the tools WILL NOT
#   be installed unless the user explicitly runs `make tools`.

ifdef ON_COCALC
  # We are on CoCalc where we install everything in the project-wide bin directory
  BIN ?= ~/.local/bin
else
  # otherwise we use a local folder (unless the user specifies BIN)
  BIN ?= .local/bin
endif

# We use sed to change the install location
MICROMAMBA_FILTER ?= sed "s:~/.local/bin:$(dir $@):g" | \
                     sed 's:YES="yes":YES="no":g'
$(BIN)/micromamba: curl
	curl -L micro.mamba.pm/install.sh | $(MICROMAMBA_FILTER) | $(SHELL)

# Everything else is installed with pipx
PIPX_PRE = PIPX_HOME=envs/ PIPX_BIN_DIR=$(dir $@)
PIPX_EXE = pipx
PIPX = $(PIPX_PRE) $(PIPX_EXE)

# The use of two targets here, the first with a requirement file, should allow make to
# use the requirements file iff it exists, falling back to the simple install.
$(BIN)/%: .tools/requirements.%.txt pipx 
	$(PIPX) install $*
	$(PIPX) inject $* $(cat $< | sed -e 's/#.*//' | tr "\n" " ")

$(BIN)/%: pipx
	$(PIPX) install $*

# In principle, we might be able to roll these into the catch-all rule,
# but it is better to be explicit here.  Without care, that will lead
# to infinite recursion.
condax:
	@command -v $@ || make $(BIN)/condax

micromamba:
	@command -v $@ || make $(BIN)/micromamba

jupytext:
	@command -v $@ || make $(BIN)/jupytext

.PHONY: condax micromamba jupytext
.NOTINTERMEDIATE: pipx   # See Notes.md.

# Note: we cannot customize the location here yet:
#    https://github.com/mariusvniekerk/condax/issues/16
# Also, there are issues on Mac OS X:
#    https://github.com/mariusvniekerk/condax/issues/63
~/.local/bin/anaconda-project: condax
	condax install anaconda-project
	condax inject anaconda-project anaconda-client

# Last-Resort default rule to check if a command exists.
%::
	@command -v $@ || \
  echo "I do not know how to make \"$@\". If it is a command, please install (e.g. with apt-get)."

$(ENV_TOOLS): .tools/environment.tools.yaml
	$(_MICROMAMBA) create -yp $@ -f $<

# Here we provide some fallback targets, but these are really not very good.
# Tools should be installed locally where possible.
USER_INSTALL_OK ?= false
ifeq ($(USER_INSTALL_OK), true)

pipx: python3
	@command -v $@ || pip install --user pipx

.PHONY pipx

endif


ifdef BIN
  export PATH := $(BIN):$(PATH)
endif

~/.local/bin/mmf_setup:
	python3 -m pip install --user --upgrade mmf-setup
ifdef COCALC_ANACONDA
	mmf_setup cocalc
endif # cookiecutter.make_tools
# ------- Documentation  -----
DOC_REQUIREMENTS =

html: $(ENV_PATH) $(DOC_REQUIREMENTS)
	$(_RUN) make -C $(DOCS) html

# We always rebuild the index.md file in case it literally includes the top-level README.md file.
# However, I do not know a good way to pass these to sphinx-autobuild yet.
ALWAYS_REBUILD ?= $(shell find $(DOCS) -type f -name "*.md" -exec grep -l '```{include}' {} + )

doc-server: $(ENV_PATH) $(DOC_REQUIREMENTS)
ifdef COCALC_ANACONDA
	$(_RUN) sphinx-autobuild --re-ignore '_build|_generated' $(DOCS) $(DOCS)/_build/html --host 0.0.0.0 --port 8000
else
	$(_RUN) sphinx-autobuild --re-ignore '_build|_generated' $(DOCS) $(DOCS)/_build/html
endif

581-Docs.tgz: $(DOCS)/*
	@make html
	tar -s "|$(DOCS)/_build/html|581-Docs|g" -zcvf $@ $(DOCS)/_build/html

581-Solutions.tgz: $(DOCS)/*
	@make html
	tar -s "|$(DOCS)/_build/html|581-Docs|g" -zcvf $@ $(DOCS)/_build/html

.PHONY: html doc-server

# ------- Experimental targets  -----
update-cookiecutter:
	# See https://github.com/cookiecutter/cookiecutter/issues/1176
	rm -rf $(DOCS)/_templates/ $(DOCS)/_static/
	$(COOKIECUTTER) --config-file $(COOKIECUTTER_YAML) --overwrite-if-exists --no-input

hg-update-cookiecutter:
	hg update cookiecutter-base
	@make update-cookiecutter
	hg commit --addremove -m "BASE: Updated cookiecutter skeleton"
	hg update default
	hg merge cookiecutter-base
	hg commit -m "Merge in cookiecutter updates"

hg-amend-cookiecutter:
	hg update cookiecutter-base
	@make update-cookiecutter
	hg amend --addremove


.PHONY: update-cookiecutter hg-update-cookiecutter hg-amend-cookiecutter

# ------- Auxilliary targets  -------
MINICONDA_SH = https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
MINICONDA_HASH = 1ea2f885b4dbc3098662845560bc64271eb17085387a70c2ba3f29fff6f8d52f

$(MINICONDA):
	wget  $(MINICONDA_SH) -qO /tmp/_miniconda.sh
	echo "$(MINICONDA_HASH)  /tmp/_miniconda.sh" > /tmp/_miniconda.shasum
	shasum -a 256 -c /tmp/_miniconda.shasum && bash /tmp/_miniconda.sh -b -p $@
	rm /tmp/_miniconda.sh*
	$@/bin/conda update -y conda
	$@/bin/conda install -y anaconda-project

	# Dropping defaults allows this to work with < 1GB
	$@/bin/conda install --override-channels --channel conda-forge -y mamba
	$@/bin/conda clean -y --all

_ext/Resources:
	-$(GET_RESOURCES)
	@if [ ! -d "$@" ]; then \
	  echo "$$RESOURCES_ERROR_MESSAGE"; \
	fi

$(DOCS)/environment.yaml: anaconda-project.yaml Makefile
	$(ANACONDA_PROJECT) run export 1> $@

# ----- Usage -----
define HELP_MESSAGE

This Makefile provides several tools to help initialize the project.  It is primarly designed
to help get a CoCalc project up an runnning, but should work on other platforms.

Variables:
   COCALC_ANACONDA: (= "$(COCALC_ANACONDA)")
                     If defined, then we assume we are on CoCalc and use this to activate
                     the conda base envrionment.  The default points to the latest version
                     of anaconda: ANACONDA_CURRENT, ANACONDA2022, ANACONDA2021, or ANACONDA2020.
                     Overrride if you want a earlier version.

                     If this is not defined, you must make sure that the ACTIVATE command
                     works properly on your machine.
   COCALC_OPTION: (= "$(COCALC_OPTION)")
                     Method for activating conda and related tools on CoCalc. Supported options:
                     * "micromamba": Install micromamba in `~/.local/bin` and use this.
                     * "anaconda": Use the version of mamba in `COCALC_ANACONDA`.
   ANACONDA_PROJECT: (= "$(ANACONDA_PROJECT)")
                     Command to run the `anaconda-project` command.  If you need to first
                     activate an environment (as on CoCalc), then this should do that.
                     Defaults to `anaconda-project`.
   ANACONDA_PROJECT_PRE: (= "$(ANACONDA_PROJECT_PRE)")
                     Pre-commands (like setting `CONDA_EXEC=mamba`) to be run before commands
                     executed with `$(ANACONDA_PROJECT)`.
   DOCS: (= "$(DOCS)")
                     Name of the documentation directory.
                     Defaults to `Docs`.
   ENV: (= "$(ENV)")
                     Name of the conda environment user by the project.
                     (Customizations have not been tested.)
                     Defaults to `phys-581`.
   ENV_PATH: (= "$(ENV_PATH)")
                     Path to the conda environment user by the project.
                     (Customizations have not been tested.)
                     Defaults to `envs/$$(ENV)`.
   ACTIVATE_ENV: (= "$(ACTIVATE_ENV)")
                     Command to activate the project environment in the shell.
                     Defaults to `conda activate -p $$(ENV_PATH)`.

Tools: By default, when not on CoCalc, we do not install tools, but have provisions here
       for doing so if desired. These variables affect the installation of these tools.

   BIN: (= "$(BIN)")
                     Location in which to install tools like `micromamba`.
   USER_INSTALL_OK: (= "$(USER_INSTALL_OK)")
                     Set to `true` if it is okay to install missing tools globally in the
                     user home directory. E.g. with `python3 -m pip install --user {}`.
                     If not `true`, then missing tools will raise an error.
   MICROMAMBA: (= "$(MICROMAMBA)")
                     Path to `micromamba`.
   SHELL_INIT_FILE: (= "$(SHELL_INIT_FILE)")
                     Init file loaded when starting the shell.  Default is `.init-file.bash`
                     but is set based on the value of `SHELL`.
   INIT_TOOLS: (= "$(INIT_TOOLS)")
                     By default, we do not install any tools with `make init`, relying on the
                     user to explicitly call `make tools`.  In some cases, however, tools
                     are always needed (CoCalc might need `micromamba` for example).  Tools
                     specified here will be installed when `make init` is called.  This can be
                     set to `tools` if all tools are desired.
   MICROMAMBA_FILTER: (= "$(MICROMAMBA_FILTER)")
                     Sed filter to install `micromamba` in the correct place and avoid
                     modifying the `.bashrc` file (run `micromamba shell init` if needed).
                     Should not need modification, but can be if the structure of the
                     install script changes.

Experimental Variables: (These features are risky or have not been full tested.)
   COOKIECUTTER_URL: (= "$(COOKIECUTTER_URL)")
                     Location of source project for cookiecutter skeletons.  Usually this is
                     `git+https://gitlab.com/forbes-group/cookiecutters.git` but can point to
                     a local directory if you have a clone or are testing changes.
   COOKIECUTTER: (= "$(COOKIECUTTER)")
                     Cookiecutter command, including `--directory` if needed.
   COOKIECUTTER_YAML: (= "$(COOKIECUTTER_YAML)")
                     Local cookiecutter yaml file for the project.

Initialization:
   make shell        Update all environments, then spawn a shell that can be used to run
                     commands in the project such as `jupyter notebook`, or `pytest`.
                     This is most similar  to `poetry shell`.  Depends on `make init`.
   make qshell       Like `make shell`, but do this quickly without checking.  This will
                     create the environment if it does not exist (`make dev`) but may not
                     perform a complete check if everything is up to date.
   make init         Initialize the environment and kernel.  On CoCalc we do specific things
                     like install mmf-setup, and activate the environment in ~/.bash_aliases.
                     This is done by `make init` if COCALC_ANACONDA is defined. If lock files
                     are provided, then the environment will be initialized from these for
                     reproducibility.
   make lock         Update the software and generate lock files by running
                     `anaconda-project lock`, `pdm lock`, etc.
   make info         Print some information.

Testing:
   make test         Runs the general tests.

Maintenance:
   make clean        Delete the documentation.
   make cleandocs    Remove documentation build.
   make cleanspace   Save disk space by removing __pycache__, calling conda clean --all etc.
   make reallyclean  Delete the documentation, environments, and kernel as well.
   make hg-update-cookiecutter 
                     Update the base branch with any pushed cookiecutter updates.  Note: this
                     assumes several things, including that you have a `default` and
                     `cookiecutter-base` base branch, as discussed in the docs, that you are
                     using mercurial, and will attempt to automatically merge the changes.
                     You may need to intervene, so try a few times manually before using this.
   make hg-amend-cookiecutter (EXPERIMENTAL)
                     Run hg amend rather than commit and does not merge
   make update-cookiecutter (EXPERIMENTAL)
                     Do the actual cookiecutter update.  Assumes appropriate VCS switching
                     and commits will be taken care of

Documentation:
   make html         Build the html documentation in `$$(DOCS)/_build/html`
   make doc-server   Build the html documentation server on http://localhost:8000
                     Uses Sphinx autobuild
   581-Docs.tgz
                     Package documentation for upload to Canvas.
endef
export HELP_MESSAGE


define RESOURCES_ERROR_MESSAGE

*************************************************************
WARNING: The `_ext/Resources` folder could not be created with

  $(GET_RESOURCES)

Likely this is because this repository is private and requires registration in the class.
If you believe that you should have access, please contact your instructor.

These resources are not crucial for the project, but are important for the course.
*************************************************************

endef
export RESOURCES_ERROR_MESSAGE
