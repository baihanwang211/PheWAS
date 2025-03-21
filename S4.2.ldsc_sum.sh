rm all.rg

FILES=($(ls | grep -E '^(scz|mdd)_.*_(EAS|EUR).*\.log$'))

for I in ${FILES[@]}; do

        # subset log files to relevant output
        tail -n5 $I | head -2 > $I.rg     # (adapt as necessary)
        
        # add to single data set
        if [[ $I == ${FILES[0]} ]]; then
            cat $I.rg > all.rg      # only including the header for the first phenotypes
        else
            cat $I.rg | sed '1d' >> all.rg
        fi

        rm $I.rg
done

