import argparse
from pipeline import MetaFluxProcessor, BiGGIDFetcher, ModuleAssigner, HMRToHumanGEMConverter
from utils import batch_fetch_bigg
'''
Guillaume PETIT IRCM Internship 2025 
Supervisor : Pr Jacques Colinge

Abstract of centent : This script as for purpose to automate the data processing of METAFlux, Compass, scFEA in order to 
                      benchmark the results of those tools

Developed by : Guillaume PETIT
Maintener : Guillaume PETIT

-----------------------------------------------------------------------

Intermediate argument are defined in order to troubleshoot the pipeline.
Please use the following command to run the script :

%bash%
python3 data_processing.py --full --input_file=<PATH_TO_METAFLUX_RAW_OUTPUT>.csv --module_file=<PATH_TO_MODULE_GENE_ASSOCIATION_TABLE_LONG_FORMAT>.csv --output_file=<PATH_TO_OUTPUT_NAME>.csv
note : if output file doesn't exist, it will be created
'''
def main():
    parser = argparse.ArgumentParser(description="Pipiline treatment for Metaflux output")
    parser.add_argument("--full", action="store_true", help="treat raw output")
    parser.add_argument("--process_metaflux", action="store_true", help="ID extraction")
    parser.add_argument("--convert_hmr", action="store_true", help="Convert HMR IDs to HumanGEM IDs")
    parser.add_argument("--fetch_bigg", action="store_true", help="fetch BiGG ID from HumanGEM IDs via API MetabolicAtlas")
    parser.add_argument("--assign_modules", action="store_true", help="match BiGG ID to scFEA metabolic modules using given association table")
    
    parser.add_argument("--input_file", type=str, required=True, help="Metaflux output file : expected format csv")
    parser.add_argument("--output_file", type=str, required=True, help="Output file for the processed dat : csv")
    parser.add_argument("--module_file", type=str, help="Association table between BiGG IDs and scFEA metabolic modules : expected format csv")
    
    args = parser.parse_args()
    
    if args.process_metaflux:
        processor = MetaFluxProcessor(args.input_file, args.output_file)
        processor.process()
    
    if args.convert_hmr:
        converter = HMRToHumanGEMConverter(args.input_file, args.output_file)
        converter.convert()
    
    if args.fetch_bigg:
        fetcher = BiGGIDFetcher(args.input_file, args.output_file)
        fetcher.fetch_bigg_ids()
    
    if args.assign_modules:
        if not args.module_file:
            parser.error("--assign_modules nécessite --module_file")
        assigner = ModuleAssigner(args.input_file, args.module_file, args.output_file)
        assigner.assign_modules()
    
    if args.full:
        intermediate_file_1 = "intermediate_metaflux.csv"  # Metaflux transform
        intermediate_file_2 = "intermediate_bigg_mapping.csv"  # BiGG ID fetcher

        # step 1 : process Metaflux output
        processor = MetaFluxProcessor(args.input_file, intermediate_file_1)
        processor.process()

        # step 2 : fetch BiGG IDs
        batch_fetch_bigg(intermediate_file_1, intermediate_file_2)

        # step 3 : assign modules
        assigner = ModuleAssigner(intermediate_file_2, args.module_file, args.output_file)
        assigner.assign_modules()

        print("Done")
    
if __name__ == "__main__":
    main()