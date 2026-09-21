rule train_known_microexons_combine_and_split_data:
  input:
    known_foreground_events=expand(
      "resources/inputs/novel_microexons/MicEvents.blacklist_filtered.{species}.tab.gz",
      species=config['species']
    ),
    known_background_events=expand(
      "resources/inputs/novel_microexons/"
      "cryptic_microexons.match_distribution.blacklist_filtered.{species}.tab.gz",
      species=config['species']
    ),
    test_set_event_ids="resources/inputs/novel_microexons/test_set_event_ids_multispecies.txt"
  output:
    training_set_foreground="resources/inputs/train_model/training_set_foreground.tab.gz",
    training_set_background="resources/inputs/train_model/training_set_background.tab.gz",
    test_set_foreground="resources/inputs/train_model/test_set_foreground.tab.gz",
    test_set_background="resources/inputs/train_model/test_set_background.tab.gz",
    validation_set_foreground="resources/inputs/train_model/validation_set_foreground.tab.gz",
    validation_set_background="resources/inputs/train_model/validation_set_background.tab.gz",
  params:
    validation_split=0.1
  conda: "../envs/base_env.yml"
  script: "../scripts/train_model/known_microexons_combine_and_split_data.py"


rule train_known_microexons_combine_and_split_data_human_only:
  input:
    known_foreground_events="resources/inputs/novel_microexons/MicEvents.blacklist_filtered.hg38.tab.gz",
    known_background_events="resources/inputs/novel_microexons/"
      "cryptic_microexons.match_distribution.blacklist_filtered.hg38.tab.gz",
    test_set_event_ids="resources/inputs/novel_microexons/test_set_event_ids_multispecies.txt"
  output:
    training_set_foreground="resources/inputs/known_microexons_hg38/training_hg38_set_foreground.tab.gz",
    training_set_background="resources/inputs/known_microexons_hg38/training_hg38_set_background.tab.gz",
    test_set_foreground="resources/inputs/known_microexons_hg38/test_hg38_set_foreground.tab.gz",
    test_set_background="resources/inputs/known_microexons_hg38/test_hg38_set_background.tab.gz",
    validation_set_foreground="resources/inputs/known_microexons_hg38/validation_hg38_set_foreground.tab.gz",
    validation_set_background="resources/inputs/known_microexons_hg38/validation_hg38_set_background.tab.gz",
  params:
    validation_split=0.1
  conda: "../envs/base_env.yml"
  script: "../scripts/train_model/known_microexons_combine_and_split_data.py"


rule train_model_known_mics:
  input:
    training_data="results/features/train_model/features_training_set.buffer",
    test_data="results/features/train_model/features_test_set.buffer",
    validation_data="results/features/train_model/features_validation_set.buffer",
    test_events_fg=rules.train_known_microexons_combine_and_split_data.output.test_set_foreground,
    test_events_bg=rules.train_known_microexons_combine_and_split_data.output.test_set_background
  output:
    test_predictions_out="resources/models/known_microexons/test_predictions.csv.gz",
    training_stats_out="resources/models/known_microexons/training_stats.json",
    xgboost_model_out="resources/models/known_microexons/xgboost_model.json"
  params:
    max_epochs=1000,
    early_stopping_rounds=100,
    xgboost_params=dict(
      eta=0.1,
      max_depth=4,
      max_delta_step=7,
      gamma=.5
    ),
    scale_pos_weight=True
  conda: "../envs/xgboost.yml"
  script: "../scripts/train_model/train_model.py"


rule train_model_known_mics_human_only:
  input:
    training_data="results/features/known_microexons_hg38/features_training_hg38_set.buffer",
    test_data="results/features/known_microexons_hg38/features_test_hg38_set.buffer",
    validation_data="results/features/known_microexons_hg38/features_validation_hg38_set.buffer",
    test_events_fg=rules.train_known_microexons_combine_and_split_data_human_only.output.test_set_foreground,
    test_events_bg=rules.train_known_microexons_combine_and_split_data_human_only.output.test_set_background
  output:
    test_predictions_out="resources/models/known_microexons_hg38/test_predictions.csv.gz",
    training_stats_out="resources/models/known_microexons_hg38/training_stats.json",
    xgboost_model_out="resources/models/known_microexons_hg38/xgboost_model.json"
  params:
    max_epochs=1000,
    early_stopping_rounds=100,
    xgboost_params=dict(
      eta=0.1,
      max_depth=4,
      max_delta_step=7,
      gamma=.5
    ),
    scale_pos_weight=True
  conda: "../envs/xgboost.yml"
  script: "../scripts/train_model/train_model.py"

rule novel_microexons_combine_and_split:
  input:
    known_foreground_events=expand(
      "resources/inputs/novel_microexons/MicEvents.blacklist_filtered.{species}.tab.gz",
      species=config['species']
    ),
    known_background_events=expand(
      "resources/inputs/novel_microexons/"
      "cryptic_microexons.match_distribution.blacklist_filtered.{species}.tab.gz",
      species=config['species']
    ),
    novel_microexons="resources/novel_microexon_hg38/hannes_total.2M.hg38.neuronal_expanded.txt",
    test_set_event_ids="resources/inputs/novel_microexons/test_set_event_ids_multispecies.txt"
  output:
    training_set_foreground="resources/inputs/novel_microexons/training_set_foreground.tab.gz",
    training_set_background="resources/inputs/novel_microexons/training_set_background.tab.gz",
    test_set_foreground="resources/inputs/novel_microexons/test_set_foreground.tab.gz",
    test_set_background="resources/inputs/novel_microexons/test_set_background.tab.gz",
    validation_set_foreground="resources/inputs/novel_microexons/validation_set_foreground.tab.gz",
    validation_set_background="resources/inputs/novel_microexons/validation_set_background.tab.gz",
  params:
    validation_split=0.1,
    test_chroms=["chr1", "chr3", "chr5", "chr7", "chr9"]
  conda: "../envs/base_env.yml"
  script: "../scripts/train_model/novel_mics_combine_and_split_data.py"

rule train_prepare_data:
  input:
    events_fg="resources/inputs/{model}/{set}_set_foreground.tab.gz",
    events_bg="resources/inputs/{model}/{set}_set_background.tab.gz",
    features_fg=expand(
      "results/features/novel_microexons/"
      "features.MicEvents.blacklist_filtered.{species}.{feature}.csv.gz",
      species=config['species'],
      feature=config['model_features']
    ) +
    expand(
      "results/features/novel_microexons/features.novel_neuronal_microexons.expanded.hg38.{feature}.csv.gz",
      feature=config['model_features']
    ),
    features_bg=expand(
      "results/features/novel_microexons/"
      "features.cryptic_microexons.match_distribution.blacklist_filtered.{species}.{feature}.csv.gz",
      species=config['species'], feature=config['model_features']
    )
  params:
    feature_names=config["model_features"]
  output:
    temp("results/features/{model}/features_{set}_set.buffer")
  conda: "../envs/xgboost.yml"
  script: "../scripts/train_model/prepare_data.py"


rule novel_mics_train_model:
  input:
    training_data="results/features/novel_microexons/features_training_set.buffer",
    test_data="results/features/novel_microexons/features_test_set.buffer",
    validation_data="results/features/novel_microexons/features_validation_set.buffer",
    test_events_fg=rules.novel_microexons_combine_and_split.output.test_set_foreground,
    test_events_bg=rules.novel_microexons_combine_and_split.output.test_set_background
  output:
    test_predictions_out="resources/models/novel_microexons/test_predictions.csv.gz",
    training_stats_out="resources/models/novel_microexons/training_stats.json",
    xgboost_model_out="resources/models/novel_microexons/xgboost_model.json"
  params:
    max_epochs=1000,
    early_stopping_rounds=100,
    xgboost_params=dict(
      eta=0.1,
      max_depth=4,
      max_delta_step=7,
      gamma=.5
    ),
    scale_pos_weight=True
  conda: "../envs/xgboost.yml"
  script: "../scripts/train_model/train_model.py"
