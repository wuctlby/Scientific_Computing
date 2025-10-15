# Scientific_Computing

This is a repository for the course of Scientific Computing

---

## 01: Docker

Brief introduction for quickly using the Docker tools for Windows

- WSL
- Docker
- VS Code

---

## 02: Computational hello word

Implement computational "hello word" meme with interepted and compiled language

### Intereted language - Python

- Create an python enviroment and install dependencies by `uv`

    Step 1: create an enviroment
    ```bash
    python3 -m venv .venv # `.venv` must be used to satisfy the default env name of `uv`
    source .venv/bin/activate
    ```

    Step 2: install uv and needed dependencies
    ```bash
    pip install uv # install uv, the modden project manager tool
    uv init # initialize the project --> produce a `pyproject.toml`
    uv add numpy # or add 'numpy'  into [project.dependencies] part of `pyproject.toml`
    python -m pip show numpy # confirm we install the `numpy` successfully
    ```
- Run `Task02_Computational_Hello_World/computational_hello_world.py` to perform the computational "hello word"
    ```bash
    # compare the performence of standard `Python List` and `Numpy` arrays
    python3 ./Task02_Computational_Hello_World/computational_hello_world.py
    ```





