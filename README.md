# PhIP_nextflow

### Docker

The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

Follow the steps below to setup your docker image:

*Note: [docker](https://www.docker.com/) is a prerequisite.*

```bash
docker build -t finterly/phip-nextflow .
```

And you're done! To run the pipeline, simply add `-profile docker`. 

```bash
nextflow run processPhage.nf -profile docker
```

