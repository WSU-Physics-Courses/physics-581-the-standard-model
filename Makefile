# Modeled somewhat after
# https://github.com/simoninireland/introduction-to-epidemics/blob/master/Makefile
SHELL = /bin/bash
_SHELL = $(notdir $(SHELL))
DOCS ?= Docs
COOKIECUTTER_URL ?= git+https://gitlab.com/forbes-group/cookiecutters.git

COOKIECUTTER ?= cookiecutter $(COOKIECUTTER_URL) --directory project
COOKIECUTTER_YAML ?= .cookiecutter.yaml

GET_RESOURCES = git clone git@gitlab.com:wsu-courses/physics-581-the-standard-model_resources.git _ext/Resources

ENV ?= phys-581
ENV_TOOLS ?= tools
ENVS ?= envs
ENV_PATH ?= $(abspath $(ENVS)/$(ENV))
ENV_TOOLS_PATH ?= $(abspath $(ENVS)/$(ENV_TOOLS))

INSTALL_TOOLS_OK ?= true
LOCAL_INSTALLER ?= micromamba

# Customize extras here pip install .[$(EXTRAS)]
EXTRAS ?= all

# Set if you want editable installs for developing
EDITABLE ?= false

BIN = $(LOCAL)/bin
MAMBA_ROOT_PREFIX ?= $(LOCAL)/micromamba
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
	# * On CoCalc, the default is to use BIN=~/.local/bin so that tools can be shared.
  #   These WILL be installed with `make init` etc.
  INIT_DEPS += $(BIN)/mmf_setup
	USER_INSTALL_OK ?= true
  ifeq ($(COCALC_OPTION), anaconda)
    # Old approach using anaconda-project in the ANACONDA_CURRENT environment.
    # Due to the $(ANACONDA_CURRENT)/.condarc issue, we must use mamba in this case
    # https://github.com/Anaconda-Platform/anaconda-project/issues/334#issuecomment-911918761
    MAMBA_EXE ?= $(ANACONDA_CURRENT)/bin/mamba
    CREATE_ENV ?= $(CONDA_PRE) $(CONDA_EXE) env create -y -f $< -p $@
    UPDATE_ENV ?= $(CONDA_PRE) $(CONDA_EXE) env update -y -f $< -p $@
    ACTIVATE_ENV ?= source $(ANACONDA_CURRENT)/bin/activate && $(CONDA_EXEC) activate -p $(ENV_PATH)
    ANACONDA_PROJECT_PRE ?= CONDA_EXE=$(CONDA_EXE) CONDA_ROOT=$(ANACONDA_CURRENT)
  else
    # New approach - use our own installation of micromamba
    CONDA_PRE ?= CONDA_EXE=$(MAMBA_EXE) CONDA_ROOT=$(LOCAL)

    #export MAMBA_ROOT_PREFIX = ~/micromamba
    #MAMBA_PRE = MAMBA_ROOT_PREFIX=$(MAMBA_ROOT_PREFIX)
    UPDATE_ENV ?= $(_MICROMAMBA) update -y -f $< -p $@
    ACTIVATE_ENV ?= $(_MICROMAMBA) activate $(ENV_PATH)
    INIT_TOOLS += $(MAMBA_EXE)
  endif
  CONDA_EXE ?= $(MAMBA_EXE)
else
	USER_INSTALL_OK ?= false
endif

ifeq ($(USER_INSTALL_OK), true)
  # Local install allowed
  LOCAL ?= ~/.local
  PYTHON3 ?= python3
  PIP_ARGS = --user --no-warn-script-location
else
  # otherwise we use a local folder (unless the user specifies BIN)
  LOCAL ?= .local
  PYTHON3 ?= $(BIN)/python3
endif

######## BIN MUST BE DEFINED BY HERE

######################################################################
# Micromamba
# We use micromamba to install python if needed.  It may also be used
# to manage environments.

# This is to resolve Issue #1.
CONDA_PRE += CUSTOM_CONDA_PREFIX="$(ENV_PATH)"

CONDA_PRE += MAMBA_ROOT_PREFIX=$(MAMBA_ROOT_PREFIX)

MICROMAMBA_EXE ?= micromamba
MAMBA_EXE ?= $(MICROMAMBA_EXE)
_MICROMAMBA = eval "$$($(MAMBA_EXE) shell hook --shell=$(_SHELL) 2> /dev/null)" && \
              $(CONDA_PRE) micromamba
# For some reason, we need to execute the function micromamba not the command here on
# CoCalc.  Previously micromamba was $(MAMBA_EXE) here.

MICROMAMBA_CREATE_ENV = $(_MICROMAMBA) create -y -f $< -p $@

######################################################################
# Generic commands: These are used in rules below, but customized
# for different installation methods like conda, micromamba, etc.
#
# RUN: run a command like sphinx-build in the appropriate environment
CREATE_ENV ?= $(MICROMAMBA_CREATE_ENV)

# The following allows micromamba to be run, even if the shell is not initialized.
_RUN = $(_MICROMAMBA) run -p $(ENV_PATH)

CONDA_EXE ?= conda
CONDA ?= eval "$$($(CONDA_PRE) $(CONDA_EXE) shell.$(_SHELL) hook)" && \
         $(CONDA_PRE) $(CONDA_EXE)
ACTIVATE_ENV ?= $(CONDA) activate -p $(ENV_PATH)


INIT_DEPS += _ext/Resources

ifeq ($(COCALC_OPTION), micromamba)
INIT_DEPS += $(BIN)/micromamba
endif
INIT_DEPS += environment.yaml
INIT_DEPS += pyproject.toml
INIT_DEPS += $(ENV_PATH)
# We have some build issues on the ARM environment when using Manim for example, so we use
# Rosetta to emulate the osx-64 platform.  Set USE_ARM = true to override this behaviour.
USE_ARM ?=

# Manim is not ready for ARM yet, so we use an intel build and rosetta
ifdef USE_ARM
else
 ifeq ($(shell uname -p),arm)
  CONDA_SUBDIR = osx-64
 endif
endif

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
	# BIN=$(BIN)
	# INIT_TOOLS=$(INIT_TOOLS)
	# INIT_DEPS=$(INIT_DEPS)
	# MICROMAMBA_EXE=$(MICROMAMBA_EXE)
	# MAMBA_EXE=$(MAMBA_EXE)
	$(MICROMAMBA) info $(ENV_PATH)

.PHONY: help usage info

SHELL_INIT_FILE ?= .init-file.$(_SHELL)
shell: $(ENV_PATH)
	$(_RUN) $(SHELL) --init-file $(SHELL_INIT_FILE)

.PHONY: shell


init: $(INIT_TOOLS) $(INIT_DEPS)
ifdef ON_COCALC
	if ! grep -Fq '$(ACTIVATE_ENV)' ~/.bash_aliases; then \
	  echo '$(ACTIVATE_ENV)' >> ~/.bash_aliases; \
	fi
	make sync
endif
	$(_RUN) python3 -m ipykernel install --user --name "phys-581" \
                                       --display-name "Python 3 (phys-581)"
ifeq ($(EDITABLE), true)
  PIP_INSTALL_ARGS += -e
endif

$(ENV_PATH): environment.yaml pyproject.toml
	$(CREATE_ENV)
	$(_RUN) python3 -m pip install $(PIP_INSTALL_ARGS) .[$(strip  $(EXTRAS))]

# Jupytext
JUPYTEXT_FILTER ?= 2>&1 | grep -v 'Reading\|Loading\|is not a paired notebook'
sync: jupytext
	find . -name ".ipynb_checkpoints" -prune -o \
	       -name "_ext" -prune -o \
	       -name "envs" -prune -o \
	       -name "*.ipynb" -o -name "*.md" \
	       -exec $(JUPYTEXT) --sync {} + $(JUPYTEXT_FILTER)



cleanspace:
	-find . -name "__pycache__" -exec $(RM) -r {} +
	-$(RM) -r _htmlcov .coverage .pytest_cache build
	-$(_MICROMAMBA) clean --all --yes

cleandocs:
	-$(RM) -r $(DOCS)/_build

clean: cleanspace cleandocs

reallyclean: realclean
realclean: clean
	$(RM) -r $(ENV_PATH)
	$(RM) -r envs
ifneq ($(USER_INSTALL_OK), true)
	$(RM) -r $(LOCAL)
endif

# https://stackoverflow.com/a/21188136

.PHONY: init lock sync clean cleanspace reallyclean realclean tools

test: $(ENV_PATH)
	$(_RUN) pytest

.PHONY: test

ifdef BIN
  export PATH := $(BIN):$(PATH)
endif

# We use sed to change the install location
MICROMAMBA_FILTER ?= sed "s:BIN_FOLDER=.*:BIN_FOLDER=$(dir $@):g" | \
                     sed "s:PREFIX_LOCATION=.*:PREFIX_LOCATION=$(LOCAL)/micromamba:g" | \
                     sed 's:YES="yes":YES="no":g'
$(BIN)/micromamba: curl
	curl -L micro.mamba.pm/install.sh | $(MICROMAMBA_FILTER) | $(SHELL)

# Everything else is installed with pipx
PIPX_HOME ?= envs/pipx
PIPX_PRE = PIPX_HOME=$(PIPX_HOME) PIPX_BIN_DIR=$(dir $@)
PIPX = $(PIPX_PRE) pipx

$(BIN)/pdm:
	$(PIPX) install pdm[all]

# The use of two targets here, the first with a requirement file, should allow make to
# use the requirements file iff it exists, falling back to the simple install.
$(BIN)/%: .tools/requirements.%.txt pipx
	$(PIPX) install $*
	$(PIPX) inject $* $(cat $< | sed -e 's/#.*//' | tr "\n" " ")

$(BIN)/%: pipx
	$(PIPX) install $*

JUPYTEXT ?= jupytext

# In principle, we might be able to roll these into the catch-all rule,
# but it is better to be explicit here.  Without care, that will lead
# to infinite recursion.
ifeq ($(INSTALL_TOOLS_OK), true)
#################################################

python3:
	@command -v $@ || make $(BIN)/python3

condax:
	@command -v $@ || make $(BIN)/condax

micromamba:
	@command -v $@ || make $(BIN)/micromamba

mmf_setup:
	@command -v $@ || make $(BIN)/mmf_setup

jupytext:
	@command -v $@ || make $(BIN)/jupytext

sphobjinv:
	@command -v $@ || make $(BIN)/sphobjinv

sphinx-autobuild:
	@command -v $@ || make $(BIN)/sphinx-autobuild

pdm:
	@command -v $@ || make $(BIN)/pdm

pipx: $(PYTHON3)
	@command -v $@ || ( $(PYTHON3) -m ensurepip --upgrade && \
	                    $(PYTHON3) -m pip install --upgrade $(PIP_ARGS) pip && \
	                    $(PYTHON3) -m pip install $(PIP_ARGS) $@ )

#################################################
endif

.PHONY: python3 condax micromamba jupytext
.NOTINTERMEDIATE: pipx pip   # See Notes.md.

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


# If .tools/environment.tools.yaml, we use it to build the $(LOCAL) environment,
# otherwise we use a venv.  See
# https://stackoverflow.com/a/1077707/1088938

tools: $(LOCAL)

$(LOCAL) $(BIN)/python3: $(ENV_TOOLS_PATH)
	ln -fs $(ENV_TOOLS_PATH)/bin/python3 $(BIN)/python3
	@echo 'Installed tools in $(LOCAL).  Evaluate the following to use:'
	@echo
	@echo "  export PATH=\"$(abspath $(LOCAL)/bin):"'$${PATH}"'
	@echo

ifneq ($(wildcard .tools/environment.tools.yaml),) 

$(ENV_TOOLS_PATH): .tools/environment.tools.yaml micromamba
	$(_MICROMAMBA) create -y -f $< -p $@

else ifeq ($(LOCAL_INSTALLER), micromamba)

$(ENV_TOOLS_PATH): micromamba
	$(_MICROMAMBA) create -c conda-forge -y -p $@ "python=3"

else

$(ENV_TOOLS_PATH):
	@if command -v python3; then \
	  python3 -m venv $@ && \
	  $(BIN)/python3 -m ensurepip --upgrade && \
	  $(BIN)/python3 -m pip install --upgrade $(PIP_ARGS) pip; \
	else
	  echo "No python3 to build $(ENV_TOOLS_PATH). Install or set LOCAL_INSTALLER=micromamba"; \
	  exit; \
	fi

endif 

$(BIN)/mmf_setup:
	$(PIPX) install mmf-setup
ifdef ON_COCALC
	mmf_setup cocalc
endif

$(BIN)/sphobjinv:
	$(PIPX) install sphobjinv

# ------- Documentation  -----
INIT_TOOLS += sphobjinv
DOC_REQUIREMENTS =

html: $(ENV_PATH) $(DOC_REQUIREMENTS)
	$(_RUN) make -C $(DOCS) html

# We always rebuild the index.md file in case it literally includes the top-level README.md file.
# However, I do not know a good way to pass these to sphinx-autobuild yet.
ALWAYS_REBUILD ?= $(shell find $(DOCS) -type f -name "*.md" -exec grep -l '```{include}' {} + )

doc-server: $(ENV_PATH) $(DOC_REQUIREMENTS) sphinx-autobuild
ifdef ON_COCALC
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

# This one is a little complicated in an attempt to solve issue #1.
HG_UPDATE_ALLOW_DIRTY_REPO ?= false
HG_UPDATE_BASE_MSG ?= "BASE: Updated cookiecutter skeleton"
HG_UPDATE_MERGE_MSG ?= "Merge in cookiecutter updates"
HG_UPDATE_SIMILARITY ?= 80
HG_UPDATE_AUTOMERGE ?= true
HG_UPDATE_BRANCH ?= cookiecutter-base

ifeq ($(HG_UPDATE_ALLOW_DIRTY_REPO), true)
	ID_CHECK=echo "Allowing repo with uncommited changes..." && id="$${id%+}"
else
	ID_CHECK=echo "Repo has uncommited changes: Aborting" && exit -1
endif

hg-init:
	if [[ ./ -ef $$(hg root 2> /dev/null) ]]; then echo "Found .hg! Aborting"; exit 0; fi
	hg init
	hg commit --addremove -m "Initial autocommit of project (make hg-init)"
	hg up null
	hg branch $(HG_UPDATE_BRANCH)
	hg revert -r default $(COOKIECUTTER_YAML)
	hg add $(COOKIECUTTER_YAML)
	hg commit $(COOKIECUTTER_YAML) -m "BASE: Initial $(COOKIECUTTER_YAML) (make hg-init)"
	hg up default
	@make hg-update-cookiecutter

# See https://stackoverflow.com/a/29085760/1088938
hg-update-cookiecutter:
	@{ \
		set -o xtrace ;\
		id=$$(hg id -i) ;\
		[[ "$${id}" = *+ ]] && $(ID_CHECK) ;\
		branch=$$(hg branch) ;\
		tmpdir=$$(mktemp -d -p . -t hg_tmp) ;\
		ln -s ../.hg "$${tmpdir}/.hg" ;\
		hg debugsetparents -R "$${tmpdir}" $(HG_UPDATE_BRANCH) ;\
		hg debugrebuildstate -R "$${tmpdir}" ;\
		hg branch -R "$${tmpdir}" $(HG_UPDATE_BRANCH) ;\
		cp $(COOKIECUTTER_YAML) "$${tmpdir}/$(COOKIECUTTER_YAML)" ;\
		$(COOKIECUTTER) --config-file "$${tmpdir}/$(COOKIECUTTER_YAML)" \
		                --no-input --output-dir "$${tmpdir}/" ;\
		hg addremove -R "$${tmpdir}" --similarity $(HG_UPDATE_SIMILARITY) ;\
		hg commit -R "$${tmpdir}" -m $(HG_UPDATE_BASE_MSG) ;\
		hg debugsetparents -R "$${tmpdir}" "$${id}" ;\
		hg debugrebuildstate -R "$${tmpdir}" ;\
		hg branch -R "$${tmpdir}" $${branch} ;\
		rm -r "$${tmpdir}" ;\
	}
ifeq ($(HG_UPDATE_AUTOMERGE), true)
	@echo "About to automerge... if this fails, resolve the issues and then"
	@echo "hg commit -m \"$(HG_UPDATE_MERGE_MSG)\""
	hg merge $(HG_UPDATE_BRANCH)
	hg commit -m $(HG_UPDATE_MERGE_MSG)
else
	@echo "$Don't forget to merge and comment (HG_UPDATE_AUTOMERGE=$(HG_UPDATE_AUTOMERGE)):"
	@echo "  hg merge $(HG_UPDATE_BRANCH)"
	@echo "  hg commit -m \"$(HG_UPDATE_MERGE_MSG)\""
endif

hg-amend-cookiecutter:
	hg update $(HG_UPDATE_BRANCH)
	@make update-cookiecutter
	hg amend --addremove


.PHONY: update-cookiecutter hg-update-cookiecutter hg-amend-cookiecutter hg-init

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
_ext/mathjax:
############################## REV here

$(DOCS)/environment.yaml: anaconda-project.yaml Makefile
	$(ANACONDA_PROJECT) run export 1> $@

# ----- Usage -----
define HELP_MESSAGE

This Makefile provides several tools to help initialize the project.  It is primarly designed
to help get a CoCalc project up an runnning, but should work on other platforms.

Variables:
   ON_COCALC: (= "$(ON_COCALC)")
                     If defined, then we assume we are on CoCalc.
   COCALC_OPTION: (= "$(COCALC_OPTION)")
                     Method for activating conda and related tools on CoCalc. Supported options:
                     * "micromamba": Install micromamba in `$$(BIN)=$$(LOCAL)/bin` and use this.
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
                     Name of the conda environment used by the project.
                     (Customizations have not been tested.)
                     Defaults to `phys-581`.
   ENVS: (= "$(ENVS)")
                     Folder to use for environments.  Defaults to `envs`.
   ENV_PATH: (= "$(ENV_PATH)")
                     Path to the conda environment user by the project.
                     (Customizations have not been tested.)
                     Defaults to `$$(ENVS)/$$(ENV)`.
   ACTIVATE_ENV: (= "$(ACTIVATE_ENV)")
                     Command to activate the project environment in the shell.
                     Defaults to `conda activate -p $$(ENV_PATH)`.
   EXTRAS: (= "$(EXTRAS)")
                     Extras to install with the project.  The defalt is `all`, which
                     includes everything needed to make the documentation, testing,
                     for running notebooks, etc.  If you just need to run the package,
                     Then you might like to make this `full` or `full,tests` to simplify.
   EDITABLE: (= "$(EDITABLE)")
                     Set to `true` if you want an editable install.

Tools: By default, when not on CoCalc, we do not install tools, but have provisions here
       for doing so if desired. These variables affect the installation of these tools.
       The strategy is to provide targets `$$(LOCAL)/bin/<tool>` which can be explicitly
       used to make the tools


   LOCAL: (= "$(LOCAL)")
                     Location in which to install tools like `micromamba`, `pipx` etc.
                     Binaries are installed in `$$(LOCAL)/bin`.  This is used as a `venv`
                     to install python-dependent tools unless `USER_INSTALL_OK`.
   INSTALL_TOOLS_OK: (= "$(INSTALL_TOOLS_OK)")
                     If `true`, then tools will be installed as needed in either `$$(LOCAL)/bin`
                     or the user home directory (if `USER_INSTALL_OK`).  Otherwise, the tools
                     will only be installed with `make tools`.
   USER_INSTALL_OK: (= "$(USER_INSTALL_OK)")
                     Set to `true` if it is okay to install missing tools globally in the
                     user home directory. E.g. with `python3 -m pip install --user {}`.
                     Sets `LOCAL=~/.local` unless overridden.
   INIT_TOOLS: (= "$(INIT_TOOLS)")
                     List of tools that overrides `INSTALL_TOOLS_OK = false`. Normally
                     `make init` does not install tool unless `INSTALL_TOOLS_OK = true`
                     relying on the user explicitly calling `make tools`.  In some cases,
                     however, tools are always needed (CoCalc might need `micromamba` for
                     example).  Tools specified here will be installed when `make init`
                     is called.  This can be set to `tools` if all tools are desired.
   MICROMAMBA: (= "$(MICROMAMBA)")
                     Path to `micromamba`.
   ENV_TOOLS: (= "$(ENV_TOOLS)")
                     Name of the conda environment in which to install tools and
                     python (if needed). (Customizations have not been tested.)
                     Defaults to `tools`.
   ENV_TOOLS_PATH: (= "$(ENV_TOOLS_PATH)")
                     Path to the conda environment with tools.
                     (Customizations have not been tested.)
                     Defaults to `$$(ENVS)/$$(ENV_TOOLS)`.
   SHELL_INIT_FILE: (= "$(SHELL_INIT_FILE)")
                     Init file loaded when starting the shell.  Default is `.init-file.bash`
                     but is set based on the value of `SHELL`.
   MAMBA_ROOT_PREFIX: (= "$(MAMBA_ROOT_PREFIX)")
                     Location where `micromamba` stores packages etc.
   MICROMAMBA_FILTER: (= "$(MICROMAMBA_FILTER)")
                     Sed filter to install `micromamba` in the correct place and avoid
                     modifying the `.bashrc` file (run `micromamba shell init` if needed).
                     Should not need modification, but can be if the structure of the
                     install script changes.
  PYTHON3: (= "$(PYTHON3)")
                     Can be used to overload the python interpreter.  If not provided
                     then this is either "python3" (on CoCalc) or "$$(BIN)/python3".  If
                     you want to be specific about the version, you should specify this
                     in `.tools/environment.tools.yaml`.

Experimental Variables: (These features are risky or have not been full tested.)
   COOKIECUTTER_URL: (= "$(COOKIECUTTER_URL)")
                     Location of source project for cookiecutter skeletons.  Usually this is
                     `git+https://gitlab.com/forbes-group/cookiecutters.git` but can point to
                     a local directory if you have a clone or are testing changes.
   COOKIECUTTER: (= "$(COOKIECUTTER)")
                     Cookiecutter command, including `--directory` if needed.
   COOKIECUTTER_YAML: (= "$(COOKIECUTTER_YAML)")
                     Local cookiecutter yaml file for the project.

Initialization and Tools:
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
   make realclean    Delete the documentation, environments, and kernel as well.  Note: will
                     not remove `$$(LOCAL)` if `USER_INSTALL_OK` is `true` to avoid removing
                     things outside of the project by mistake.
   make hg-update-cookiecutter 
                     Update the base branch with any pushed cookiecutter updates.
            WARNING: This should only be run on a clean repository (otherwise, uncommited files
                     might be added to `$(HG_UPDATE_BRANCH)`.)
               Note: This assumes several things, including that you have a `default` and
                     `$(HG_UPDATE_BRANCH)` base branch, as discussed in the docs, that you are
                     using mercurial, and will attempt to automatically merge the changes.
                     You may need to intervene, so try a few times manually before using this.
   make hg-amend-cookiecutter (EXPERIMENTAL)
                     Run hg amend rather than commit and does not merge
            WARNING: This should only be run on a clean repository (otherwise, uncommited files
                     might be added to `$(HG_UPDATE_BRANCH)`.)
   make update-cookiecutter (EXPERIMENTAL)
                     Do the actual cookiecutter update.  Assumes appropriate VCS switching
                     and commits will be taken care of
            WARNING: This should only be run on a clean repository (otherwise, uncommited files
                     might be added to `$(HG_UPDATE_BRANCH)`.)

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
