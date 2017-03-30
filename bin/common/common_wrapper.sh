# Verify a command line argument is not null
verify_not_null() {
   if [ -z "$1" ]; then
      usage
   fi
}

# Print usage and exit
usage()
{
    echo "Usage: $( basename $0 ) [ -site <sitename> | -params <params> ]" >&2
    exit 1
}

cleanup()
{
   if [ "$1" -eq 0 ]; then
      echo Cleaning up, please wait
      sleep 5
      rm -rf run??? finder.out logs
   fi
}
