import argparse


def main():
    parser = argparse.ArgumentParser(description="tunnel_extraction")

    # Create subparsers
    subparsers = parser.add_subparsers(dest="subcommand", title="Subcommands", description="Available subcommands", required=True)

    # Subparser for DBSCAN
    dbscan_parser = subparsers.add_parser("dbscan", help="DBSCAN-related options")
    dbscan_parser.add_argument("--eps", type=float, help="Specify eps value")
    dbscan_parser.add_argument("--min_samples", type=int, help="Specify min_samples value")
    dbscan_parser.add_argument("--metric", choices=DBSCAN_METRICS, help="Choose a metric")

    # Subparser for Plot
    plt_parser = subparsers.add_parser("plt", help="Plot-related options")
    plt_parser.add_argument("--plot", action="store_true", help="Enable plotting")
    plt_parser.add_argument("--cluster_only", type=int, help="Specify cluster only value")

    # Subparser for ETL
    etl_parser = subparsers.add_parser("etl", help="ETL-related options")
    etl_parser.add_argument("--rcsb_id", type=str, help="Specify RCSB ID")
    etl_parser.add_argument("--input", type=str, help="Specify input file path")

    args = parser.parse_args()

    # Access the values of the arguments based on the subcommand
    if args.subcommand == "dbscan":
        eps = args.eps
        min_samples = args.min_samples
        metric = args.metric
        print("DBSCAN Options:")
        print("Eps:", eps)
        print("Min Samples:", min_samples)
        print("Metric:", metric)

    elif args.subcommand == "plt":
        plot = args.plot
        cluster_only = args.cluster_only

    elif args.subcommand == "etl":
        rcsb_id = args.rcsb_id
        input_path = args.input
        print("ETL Options:")
        print("RCSB ID:", rcsb_id)
        print("Input Path:", input_path)

if __name__ == "__main__":
    main()
