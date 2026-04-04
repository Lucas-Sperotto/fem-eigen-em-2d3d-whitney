cd ..

cmake -S . -B build

cmake --build build --target helm10_rect helm10_circle helm10_coax -j2

./build/helm10_rect 1.0 14 14 20 --backend closed-form
./build/helm10_circle 10 48 20 --backend closed-form
./build/helm10_coax 10 48 20 --backend closed-form

python3 scripts/helm10.py

cmake --build build --target edge_rect edge_circle edge_coax -j2

./build/edge_rect 14 14 20 --backend closed-form
./build/edge_circle 10 48 20 --backend closed-form
./build/edge_coax 10 48 20 --backend closed-form

python3 scripts/helmvec.py


cmake --build build --target mixed_rect mixed_circle mixed_coax -j2

./build/mixed_rect 12 6 --backend closed-form
./build/mixed_circle 10 48 --backend closed-form
./build/mixed_coax 10 48 --backend closed-form

python3 scripts/helmvec1.py
