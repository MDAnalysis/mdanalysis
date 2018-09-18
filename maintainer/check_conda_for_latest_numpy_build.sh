#! /bin/bash

if [ -n "$DEBUG" ]; then
    set -x
fi

# Sometimes, conda fails due to a CondaHTTPError.
# In that case, we want to automatically retry for a certain number of times:
if [ -z "$CONDA_HTTP_ERROR_MAX_RETRY" ]; then
    CONDA_HTTP_ERROR_MAX_RETRY=3
fi

# Checking for the latest stable numpy build only makes sense if we haven't already installed it:
if [ "$NUMPY_VERSION" == "stable" ]; then
    echo "INFO: NUMPY_VERSION is set to \"stable\", skipping check for latest stable numpy build." 1>&2
    exit
fi

# Checking for the latest stable numpy build only makes sense if a specific Python version is required:
if [ -z "$PYTHON_VERSION" ]; then
    echo "INFO: PYTHON_VERSION is not set, skipping check for latest stable numpy build." 1>&2
    exit
fi

# Make sure ci-helpers actually exported LATEST_NUMPY_STABLE:
if [ -z "$LATEST_NUMPY_STABLE" ]; then
    echo "ERROR: Unable to check for latest numpy build: Environment variable LATEST_NUMPY_STABLE is not set." 1>&2
    echo "Possible reasons:"
    echo "    1. ci-helpers failed to export LATEST_NUMPY_STABLE"
    echo "    2. You executed this script before calling \"source ci-helpers/travis/setup_conda.sh\"." 1>&2
    exit 1
fi

# Make sure CONDA_CHANNELS is unset (otherwise, conda will choke on that):
if [ -n "$CONDA_CHANNELS" ]; then
    CONDA_CHANNELS=""
fi

# Make sure conda is in PATH:
command -v conda >/dev/null 2>&1
if [ "$?" != "0" ]; then
    echo "ERROR: Command \"conda\" not found. Make sure the conda bin directory is added to PATH." 1>&2
    exit 1
fi

# Get output of "conda info numpy=<version>=<build_string>", retry on CondaHTTPError
retry=0
exitval=0
conda_call="conda info numpy=$LATEST_NUMPY_STABLE=py${PYTHON_VERSION/./}*"
while true; do
    numpy_conda_info="$($conda_call 2>&1)"
    exitval="$?"
    if [ "$exitval" == "0" ]; then
        break
    fi
    if [ -n "$(echo "$numpy_conda_info" | grep 'CondaHTTPError')" ]; then
        if [[ $retry -lt $CONDA_HTTP_ERROR_MAX_RETRY ]]; then
            echo "WARNING: The comand \"$conda_call\" failed due to a CondaHTTPError, retrying."
            retry=$(($retry + 1))
            sleep 2
        else
            echo "ERROR: The following command failed due to a CondaHTTPError after $CONDA_HTTP_ERROR_MAX_RETRY retries:" 1>&2
            echo "    $conda_call" 1>&2
            echo "Output of failed command:" 1>&2
            echo "$numpy_conda_info" 1>&2
            exit 1
        fi
    else
        echo "ERROR: The following command returned with non-zero exit status $exitval:" 1>&2
        echo "    conda info numpy=$LATEST_NUMPY_STABLE=py${PYTHON_VERSION/./}*" 1>&2
        echo "Output of failed command:" 1>&2
        echo "$numpy_conda_info" 1>&2
        exit 1
    fi
done

# If the conda info output isn't empty, let's check it:
if [ -n "$numpy_conda_info" ]; then
    numpy_builds="$(echo "$numpy_conda_info" | grep "numpy *${LATEST_NUMPY_STABLE/./\\.}")"
    if [ -n "$numpy_builds" ]; then
        echo "Latest stable numpy ($LATEST_NUMPY_STABLE) build(s) available for Python $PYTHON_VERSION:"
        echo "$numpy_builds"
        echo "Please remove NUMPY_VERSION=$NUMPY_VERSION in .travis.yml for Python $PYTHON_VERSION tests! Aborting."
        exit 1
    fi
fi

# If we end up here, there's no latest stable numpy build available.
echo "INFO: No Python $PYTHON_VERSION build of latest stable numpy ($LATEST_NUMPY_STABLE) available, continuing."
