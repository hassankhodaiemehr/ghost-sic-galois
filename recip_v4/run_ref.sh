d=$1; DIG=$2
echo "d=$d; DIG=$DIG; read(\"refine_drv.gp\");" | gp -q -s 8G > ref$d.log 2>&1
