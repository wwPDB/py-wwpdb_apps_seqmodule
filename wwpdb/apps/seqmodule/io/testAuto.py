from wwpdb.apps.seqmodule.util.Autodict import Autodict


def test1():
    I = Autodict()  # noqa: E741
    I["A"]["feature"][1][1][1] = (1, 2)
    I["A"]["feature"][1][2][1] = (1, 2)
    #
    I["A"]["feature"][2][1][1] = (1, 2)
    I["A"]["feature"][2][2][1] = (1, 2)
    I["A"]["feature"][2][3][1] = (1, 2)
    print(I)
    print(I["A"]["feature"].keys())
    print(I["A"]["feature"][1].keys())
    print(I["A"]["feature"][2].keys())
    #
    I["A"]["feature"] = {}
    print(I)
    print(I["A"]["feature"].keys())
    print(I["A"]["feature"][1].keys())
    print(I["A"]["feature"][2].keys())


if __name__ == "__main__":
    test1()
