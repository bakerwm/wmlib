#!/bin/bash

# Prepare log file
# 
# save log at: ${HOME}/work/system_log/jupyter_notebook/jupyter_notebook_20201205_140329.log

function get_logfile() {
    cur_time=$(date +%H%M%S)
    year=$(date +%Y)
    mon=$(date +%m)
    wk=$(date +%W)
    day=$(date +%d)
    outdir="${HOME}/work/system_log/jupyterlab/${year}/${mon}"
    fname="jupyter_notebook_${year}${mon}${day}_${cur_time}.log"
    logfile="${outdir}/${fname}"
    mkdir -p ${outdir}
    echo $logfile
}


function install_jupyter() {
    # install jupyter 2.2.9 (3.0.x not full supported) 
    echo "Install Jupyterlab..."
    mamba install -c conda-forge -y jupyterlab=2.2.9
    mamba install -c conda-forge -y nodejs
    mamba install -c conda-forge -y qtconsole ipywidgets

    # install extensions
    echo "Install extensions..."
    echo "Install git"
    conda install -c conda-forge -y jupyterlab-git
    jupyter labextension install @jupyterlab/git
    jupyter serverextension enable --py jupyterlab_git
    
    echo "Install github"
    jupyter labextension install @jupyterlab/github
    
    echo "Install toc"
    jupyter labextension install @jupyterlab/toc

    ## failed
    # echo "Install jupyterlab_code_formatter"
    # conda install -y jupyterlab_code_formatter
    # jupyter serverextension enable --py jupyterlab_code_formatter
    # jupyter labextension install @ryantam626/jupyterlab_code_formatter

    echo "Install jupyterlab-execute-time"
    jupyter labextension install jupyterlab-execute-time

    echo "Install jupyterlab-tailwind-theme"
    jupyter labextension install jupyterlab-tailwind-theme

    echo "Build jupyter extensions"
    jupyter lab build
}




# Check if jupyter installed
# see: hash, type, command

function has_jupyter() {
    if hash jupyter 2>/dev/null 
    then
        jupyter=$(command -v jupyter)
        >&2 echo ">> [jupyterlab] jupyter_command : [ok] $jupyter"
        echo 0
    else
        >&2 echo ">> [jupyterlab] jupyter_command : [failed] jupyter not found, see https://jupyterlab.readthedocs.io/en/stable/ for installation"
        echo 1
    fi    
}


# Check if jupyter is running
# Only for single user
# $ jupyter notebook list
# Currently running servers:
# http://*:8888/ :: /home/sam/

function jupyter_port() {
    # check port
    # [[ $(has_jupyter) -gt 0 ]] && exit 1
    log=$(jupyter notebook list | grep -v 'Currently running' | cut -d":" -f 3)
    port=${log/\//}
    echo $port
}


function jupyter_url() {
    # check url
    # [[ $(has_jupyter) -gt 0 ]] && exit 1
    url=$(jupyter notebook list | grep -v 'Currently running' | cut -d" " -f 1)
    echo $url
}



# Check if jupyter is running or not
# $ jupyter notebook list
# Currently running servers:
# http://*:8888/ :: /home/sam/

function status_jupyter() {
    # [[ $(has_jupyter) -gt 0 ]] && exit 1
    url=$(jupyter_url)
    echo -n ">> [Jupyterlab] status : "
    if [[ -z $url ]] 
    then
        echo "[stopped] ... Jupyter is not running"
    else
        echo "[running] ... Jupyter running at [$url]"
        # log=$(jupyter notebook list | grep -v 'Currently running')
        # echo $log
    fi
}



# Start Jupyter notebook
#
# Check if the Jupyter is running 

function start_jupyter() {
    # [[ $(has_jupyter) -gt 0 ]] && exit 1
    echo -n ">> [Jupyterlab] start : "
    url=$(jupyter_url)
    if [[ -z $url ]] 
    then
        logfile=$(get_logfile)
        jupyter lab &> $logfile &
        sleep 5 # wait for jupyter
        url2=$(jupyter_url)
        echo "[ok] ... Jupyter running [$url2]"
        echo "Log saved : $logfile"
    else
        echo "[failed] ... Another Jupyter is running [$url]"
    fi
}




# Kill the process, bind to the port, running jupyter
#
# TARGET_PORT=8888
# kill -9 $(lsof -n -i4TCP:$TARGET_PORT | cut -f 2 -d " ") 

function stop_jupyter() {
    # [[ $(has_jupyter) -gt 0 ]] && exit 1
    echo -n ">> [Jupyterlab] stop : "
    port=$(jupyter_port)
    url=$(jupyter_url)
    if [[ -z $port ]] ; then
        echo "[ok] ... No Jupyter is running"
    else
        kill -9 $(lsof -n -i4TCP:$port | cut -f 2 -d " ")
        sleep 5 # wait for jupyter
        echo "[ok] ... Jupyter stopped [$url]"
    fi
}


# Restart Jupyter notebook
#
# Check if the Jupyter is running, re-run

function restart_jupyter() {
    echo ">> [Jupyter notebook] restart : "
    url=$(jupyter_url)
    if [[ ! -z $url ]] ; then
        stop_jupyter # stop
    fi
    start_jupyter # start
}

# Command-line arg
[[ $# -lt 1 ]] && echo "run_jupyter.sh [status|start|stop|restart|install]" && exit 1

# # check command exists
# [[ $(has_jupyter) -gt 0 ]] && exit 1

case $1 in 
    status)  status_jupyter ;;
    start)   start_jupyter ;;
    stop)    stop_jupyter ;;
    restart) restart_jupyter ;;
    install) install_jupyter ;;
    *)       echo "[$1] unknown, choose [start|stop|restart]" ;;
esac


## EOF