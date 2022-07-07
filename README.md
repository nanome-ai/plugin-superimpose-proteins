# Nanome - Superimpose Proteins Plugin

The Superimpose Proteins plugin overlays two or more proteins in 3D space for visual comparison and calculates Root Mean Square Deviation (RMSD) values as a numeric metric on structural similarity.

Protein overlay is done either by using backbone (alpha-carbons) or by using all heavy protein atoms; and the alignment can be done on the full protein or limited to a chain.

A Fixed Reference structure is chosen, and Moving Structures are translocated to be superimposed upon the Fixed Reference. RMSD values are reported between the Fixed Reference and each Moving Structure.

## Dependencies

[Docker](https://docs.docker.com/get-docker/)


## Usage

To run Superimpose in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh [run_args]
```

### Usage

These args are passed to deploy.sh, and then forwarded to the Plugin's run.py command

```sh
usage: run.py [-h] [-a HOST] [-p PORT] [-k KEY] [-n NAME] [-v] [--write-log-file WRITE_LOG_FILE]
              [--remote-logging REMOTE_LOGGING] [-r] [-i IGNORE]

Starts a Nanome Plugin.

optional arguments:
  -h, --help            show this help message and exit
  -a HOST, --host HOST  connects to NTS at the specified IP address
  -p PORT, --port PORT  connects to NTS at the specified port
  -k KEY, --keyfile KEY
                        Specifies a key file or key string to use to connect to NTS
  -n NAME, --name NAME  Name to display for this plugin in Nanome
  -v, --verbose         enable verbose mode, to display Logs.debug
  --write-log-file WRITE_LOG_FILE
                        Enable or disable writing logs to .log file
  --remote-logging REMOTE_LOGGING
                        Toggle whether or not logs should be forwarded to NTS.
  -r, --auto-reload     Restart plugin automatically if a .py or .json file in current directory changes
  -i IGNORE, --ignore IGNORE
                        To use with auto-reload. All paths matching this pattern will be ignored, use commas to specify
                        several. Supports */?/[seq]/[!seq]
```

## License

MIT
