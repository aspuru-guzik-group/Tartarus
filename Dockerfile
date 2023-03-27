FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get upgrade

SHELL ["/bin/bash", "--login", "-c"]

WORKDIR /benchmark
COPY environment.yml .

RUN conda env create -f environment.yml
RUN conda init bash
RUN echo "conda activate tartarus" >> ~/.bashrc

RUN echo "export XTBHOME=$CONDA_PREFIX" > $CONDA_PREFIX/etc/conda/activate.d/env.sh
RUN echo "source $CONDA_PREFIX/share/xtb/config_env.bash" >> $CONDA_PREFIX/etc/conda/activate.d/env.sh

RUN pip install -i https://test.pypi.org/simple/ geodesic-interpolate
RUN pip install pytest

RUN mkdir /data
COPY . .

SHELL ["conda", "run", "-n", "tartarus", "/bin/bash", "-c"]
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "tartarus", "python", "/benchmark/benchmark.py"]
