For both formats, LNK and Dict, we compare 3 different basic algorithms

#### 1. Add all non-CSC matrices, and convert the sum to CSC `AddNonCSCthenCSC`

#### 2. Convert all non-CSC matrices to CSC and add them to an empty CSC matrix `AddNonCSCditoCSC`

#### 3. Try to construct a clever way of filling the CSC matrix `be_clever_si`

We fix the matrix size to n = 1000 and the number of matrices to be added to 4.

nnz = 500

| | Time in ms | Allocs in bytes |
|-----|-----|-----|
| AddNonCSCthenCSC LNK | 0.20302596399999995 | 132080 |
|AddNonCSCditoCSC LNK | 0.3794482360000003 | 197296 |
|**be clever si     LNK** | **0.156190813** | 40544 |
|AddNonCSCthenCSC Dict | 0.4814883580000003 | 272224 |
|*AddNonCSCditoCSC Dict* | *0.4374453029999999* | 344768 |
|be clever si      Dict | 0.8084115250000005 | 638160 |

---

nnz = 1 000


| | Time in ms | Allocs in bytes |
|-----|-----|-----|
| AddNonCSCthenCSC LNK | 0.4333999360000002 | 374032 |
| AddNonCSCditoCSC LNK | 0.43911605300000023 | 301120 |
| **be clever si     LNK** | **0.25078191200000005** |  72384 |
| AddNonCSCthenCSC Dict | 1.283850133 | 722880 |
| *AddNonCSCditoCSC Dict* | *0.7667664479999997* | 559808 |
| be clever si      Dict | 1.241192799999999 |  861568 |

---

nnz = 5 000

| | Time in ms | Allocs in bytes |
|-----|-----|-----|
| AddNonCSCthenCSC LNK | 3.9810572510000006 | 987016 |
| **AddNonCSCditoCSC LNK** | **1.7841651669999987** | 1123952 |
| be clever si     LNK | 2.3840109350000027 | 325440 |
| AddNonCSCthenCSC Dict | 8.036482746 | 3090496 |
| *AddNonCSCditoCSC Dict* | *3.0687390060000026* | 2266432 |
| be clever si      Dict | 5.936740458999996 | 2607760 |

---

nnz = 20 000

| | Time in ms | Allocs in bytes |
|-----|-----|-----|
| AddNonCSCthenCSC LNK | 45.161748425000006 | 3852768 |
| **AddNonCSCditoCSC LNK** | **11.810012035999998** | 4095584 |
| be clever si     LNK | 23.044766257000017 | 1238976 |
| AddNonCSCthenCSC Dict | 35.892576227000035 | 12116608 |
| *AddNonCSCditoCSC Dict* | *21.911943386000008* | 8432896 |
| be clever si      Dict |23.884619575999963 | 7935936 |
