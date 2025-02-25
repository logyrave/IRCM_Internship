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
        intermediate_file_1 = "intermediate_metaflux.csv"  # 🔹 Étape 1 : transformation de METAFlux
        intermediate_file_2 = "intermediate_bigg_mapping.csv"  # 🔹 Étape 2 : récupération des BiGG_ID

        # 🔹 Étape 1 : Appliquer le traitement METAFlux → Sortie Metadata, Value
        processor = MetaFluxProcessor(args.input_file, intermediate_file_1)
        processor.process()

        # 🔹 Étape 2 : Associer les BiGG_ID → Doit prendre en entrée le fichier METAFlux traité
        batch_fetch_bigg(intermediate_file_1, intermediate_file_2)

        # 🔹 Étape 3 : Assigner les modules (dernière étape du pipeline)
        assigner = ModuleAssigner(intermediate_file_2, args.module_file, args.output_file)
        assigner.assign_modules()

        print("✅ Pipeline complète exécutée avec succès !")
    
if __name__ == "__main__":
    main()