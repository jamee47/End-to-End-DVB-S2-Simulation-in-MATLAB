import argparse

def test_func(num, string, is_true):
    print(num)
    print(string)
    print(is_true)

def test_func2(string2):
    print(string2)

def main():
    
    parser = argparse.ArgumentParser(description="Testing the argparse")
    parser.add_argument("--id",type= int, required = True,
                        help="Enter the ID")
    parser.add_argument('--name', type= str, required=False, default="No name",
                        help= "Enter the Name")
    parser.add_argument("--true", action="store_true",
                        help="Testing the boolean")
    parser.add_argument("--sec", type =str, default="second function")
    args = parser.parse_args()
    test_func(num=args.id,
              string=args.name,
              is_true=args.true)
    test_func2(string2=args.sec)
    
if __name__ == "__main__":
    main()