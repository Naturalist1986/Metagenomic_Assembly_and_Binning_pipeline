#!/bin/bash
# 10_create_abundance_barplots.sh - Create barplots from abundance matrices

# Run the Python script
python3 10_create_abundance_barplots.py

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "View your plots:"
    echo "========================================="
    echo ""
    echo "# List all plots"
    echo "ls -lh /sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap/abundance_plots/*.png"
    echo ""
    echo "# View a specific plot (requires X11 forwarding or download)"
    echo "display /sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap/abundance_plots/RH_abundance_barplot.png"
    echo ""
    echo "# Or download to your local machine:"
    echo "scp moshea@moriah-gw-01:/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap/abundance_plots/*.png ."
    echo ""
fi
