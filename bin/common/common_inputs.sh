BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

# Fetch a yaml parameter
get_param() {
    $BINDIR/common/param_get.py --input $params --key $1
}

area_to_long() {
    area=$1
    if [ $area = mirca ] || [ $area = iizumi ] || [ $area = spam ]; then
        echo fixed_${area}_mask
    else
        echo dynamic_${area}_mask
    fi
}
