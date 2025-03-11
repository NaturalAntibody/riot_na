# EXAMPLE HOW TO INSTALL RIOT IN DOCKER IMAGE

FROM  python:3.11.11-slim

WORKDIR /app

# Copy the current directory contents into the container at /app
COPY ./tests/test_e2e.py /app

# GCC needed for scikit-bio installation. Can be removed after they provide wheels: https://github.com/scikit-bio/scikit-bio/issues/588
RUN apt-get update && apt-get install -y gcc

RUN --mount=type=cache,target=/root/.cache pip install riot_na

# Define environment variable
ENV BASE_PATH=/app

# Test if riot is working correctly
CMD ["python", "test_e2e.py"]
