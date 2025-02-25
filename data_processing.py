import argparse
from pipeline import MetaFluxProcessor, BiGGIDFetcher, ModuleAssigner, HMRToHumanGEMConverter
from utils import batch_fetch_bigg
'''
Guillaume PETIT IRCM Internship 2025 
Supervisor : Pr Jacques Colinge
Abstract of centent : This script as for purpose to automate the data processing of METAFlux, Compass, scFEA in order to 
                      benchmark the results of those tools
'''
def main():
    parser = argparse.ArgumentParser(description="Pipeline de traitement des données MetaFlux")
    parser.add_argument("--full", action="store_true", help="process all steps to treat raw output")
    parser.add_argument("--process_metaflux", action="store_true", help="Extraire les IDs MAR0xxxx depuis MetaFlux output")
    parser.add_argument("--convert_hmr", action="store_true", help="Convertir les IDs HMR en Human-GEM")
    parser.add_argument("--fetch_bigg", action="store_true", help="Récupérer les identifiants BiGG depuis MAR0xxxx via API MetabolicAtlas")
    parser.add_argument("--assign_modules", action="store_true", help="Associer les BiGG IDs aux modules scFEA")
    
    parser.add_argument("--input_file", type=str, required=True, help="Fichier d'entrée pour l'étape choisie")
    parser.add_argument("--output_file", type=str, required=True, help="Fichier de sortie pour l'étape choisie")
    parser.add_argument("--module_file", type=str, help="Fichier contenant la correspondance Gene_ID -> Module")
    
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
        intermediate_file = "intermediate_bigg_mapping.csv"
        
        batch_fetch_bigg(args.input_file, intermediate_file)
        
        assigner = ModuleAssigner(intermediate_file, args.module_file, args.output_file)
        assigner.assign_modules()
        print("Pipeline complète exécutée avec succès.")    
    
if __name__ == "__main__":
    main()