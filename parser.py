import json
import argparse
from pathlib import Path

def load_config(path):
    with open(path, "r") as j:
        config = json.load(j)
    return config


def get_args():
    config = load_config("config.json")
    # parent_path = Path(__name__).parent
    #
    # relative_data_path = Path(config["data_path"])
    # relative_output_path = Path(config["output_path"])
    # relative_model_path = Path(config["model_path"])

    parser = argparse.ArgumentParser(description="scrna conduction")

    # Path and file configs
    # parser.add_argument(
    #     "--data_path", default=relative_data_path, help="The dataset path.", type=str
    # )
    #
    # parser.add_argument(
    #     "--output_path",
    #     default=relative_output_path,
    #     help="The predictions will be saved to this path.",
    #     type=str,
    # )
    parser.add_argument(
        "--clustering_method", default=config["train_file"], help="default=train.tsv", type=list
    )
    parser.add_argument(
        "--val_file", default=None, help="validation file, default=None", type=str
    )
    parser.add_argument(
        "--test_file", default=config["test_file"], help="default=test.tsv", type=str
    )
    parser.add_argument(
        "--train_column",
        default=config["train_column"],
        help="The train column, default=Phrase",
        type=str,
    )
    parser.add_argument(
        "--label_column",
        default=config["label_column"],
        help="The label column, default=labels",
        type=str,
    )
    parser.add_argument(
        "--identifier_column",
        default=config["identifier_column"],
        help="The identifier column used for prediction file, default=PhraseId",
        type=str,
    )

    # Model configs
    parser.add_argument(
        "--tokenizer",
        default=config["tokenizer"],
        help="pre-trained tokenizer, default=distilbert-base-uncased",
        type=str,
    )

    parser.add_argument(
        "--embedding_dim",
        default=config["embedding_dim"],
        help="Tokens will be embedded to a vector, default=128",
        type=int,
    )
    parser.add_argument(
        "--hidden_dim",
        default=config["hidden_dim"],
        help="The hidden state dim of BiLSTM. default=256",
        type=int,
    )
    parser.add_argument(
        "--output_dim",
        default=config["output_dim"],
        help="output dim of BiLSTM-> num of label class, default=5",
        type=int,
    )
    parser.add_argument(
        "--num_layers",
        default=config["num_layers"],
        help="number of LSTM layers, default=2",
        type=int,
    )

    # Optimizer config
    parser.add_argument(
        "--lr",
        default=config["lr"],
        help="Learning rate of the optimizer,default=1e-4",
        type=float,
    )

    # Training configs
    parser.add_argument(
        "--train_batch_size",
        default=config["train_batch_size"],
        help="training batch size, default=64",
        type=int,
    )
    parser.add_argument(
        "--test_batch_size",
        default=config["test_batch_size"],
        help="test batch size, default=64",
        type=int,
    )
    parser.add_argument(
        "--num_epochs",
        default=config["num_epochs"],
        help="training epochs number, default=10.",
        type=int,
    )
    parser.add_argument(
        "--save_every_n_epoch",
        default=config["save_every_n_epoch"],
        help="checkpoint saved every n epoches, default=5",
        type=int,
    )
    parser.add_argument(
        "--world_size",
        default=config["world_size"],
        help="number of processes started, default=all gpu",
        type=int,
    )

    # Mode config
    parser.add_argument(
        "--test",
        default=config["test"],
        help="Test on the testset.",
        action="store_true",
    )
    parser.add_argument(
        "--run_from_ckp",
        default=config["run_from_ckp"],
        help="rerun from checkpoint",
        action="store_true",
    )

    args = parser.parse_args()
    args.tokenizer = AutoTokenizer.from_pretrained(args.tokenizer)
    args.vocab_size = args.tokenizer.vocab_size

    return args