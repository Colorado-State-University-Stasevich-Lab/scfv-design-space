#!/bin/bash
# Clone ProteinMPNN if not already present
if [ ! -d "ProteinMPNN" ]; then
    echo "Cloning ProteinMPNN..."
    git clone https://github.com/dauparas/ProteinMPNN.git
fi
