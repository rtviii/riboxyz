from ribctl.lib.npet2.run import run_npet2

def main():
    ctx = run_npet2("7K00")
    print("run_dir:", ctx.store.run_dir)

if __name__ == "__main__":
    main()
