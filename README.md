1. state_arrow.c 파일에 관하여

필요한 자료: DATA_truth 폴더와 DATA_matrix 폴더의 자료

컴파일: gcc -std=c99 -Wall -g -lm state_arrow.c -o state_arrow.out

실행: ./state_arrow.out {origin_filename}

예시: ./state_arrow.out Signal_Transduction_in_Fibroblasts

출력결과: ./STATE_Networks/Signal_Transduction_in_Fibroblasts.nw 파일이 생성

생성된 파일은 2열로 된 파일이며 {source} {target} 순으로 저장되었다. 방향성이 있는 네트워크에서 노드 사이의 연결을 나타냄.
노드를 의미하는 숫자는 2진법의 Boolean state를 10진법으로 환산한 값이다. 

Boolean state는 ./DATA_matrix/Signal_Transduction_in_Fibroblasts_MATRIX_pm.csv 의 첫째줄에 나와있는 유전자 순으로 0/1 값을 가지는지를 의미한다.

