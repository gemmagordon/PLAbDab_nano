import torch
import pandas as pd
import re
import numpy as np
from PLAbDab_nano.database_generate.shark_anarcii.model import Transformer


start_token, end_token, pad_token = '<SOS>', '<EOS>', '<PAD>'

class Tokenizer:
    @staticmethod
    def seq_vocab():
        aa = 'acdefghiklmnopqrstvwyX'
        aa = [x.upper() for x in aa]
        vocab = [pad_token] + [start_token] + [end_token] + aa
        return vocab

    @staticmethod
    def seq_char_to_int():
        vocab = Tokenizer.seq_vocab()   
        character_to_integer = {char: i for i, char in enumerate(vocab)}
        return character_to_integer

    @staticmethod
    def num_vocab():
        nums = [str(x) for x in range(0, 150)]
        nums.append("X")  # X for insert
        vocab = [pad_token] + [start_token] + [end_token] + nums
        return vocab

    @staticmethod
    def number_char_to_int():
        vocab = Tokenizer.num_vocab()
        character_to_integer = {char: i for i, char in enumerate(vocab)}
        return character_to_integer

    def seq_tokeniser(ls):
        character_to_integer = Tokenizer.seq_char_to_int()
        integer_encoded = np.array([character_to_integer[char] for char in ls], dtype=np.int16)
        return integer_encoded

    def number_tokeniser(ls):
        character_to_integer = Tokenizer.number_char_to_int()
        integer_encoded = np.array([character_to_integer[char] for char in ls], dtype=np.int16)
        return integer_encoded

def tokenise_sharks(ls):
    seq_ls = [[x for x in seq] for seq in ls]
    seq_ls = [["<SOS>"] + seq + ["<EOS>"] for seq in seq_ls]
    X = [Tokenizer.seq_tokeniser(seq) for seq in seq_ls]
    X = [torch.from_numpy(i).type(torch.int32).unsqueeze(0) for i in X]
    return X

def build_inward_list(length, start_num, end_num):
    result = []
    midpoint = length // 2  # Find the middle index by floor division
    # 5 becomes 2 ensures that we start at the higher number if uneven
    for i in range(midpoint):
        result.append((str(start_num), chr(ord('A') + i)))
    for i in range(midpoint, length):
        if length % 2 != 0: # odd
            result.append((str(end_num), chr(ord('A') + (midpoint*2-i))))
        elif length % 2 == 0: # even
            result.append((str(end_num), chr(ord('A') + (midpoint*2-(i+1)))))
    return result

def number(ls, device=torch.device('cpu')):
    INPUT_DIM = 25
    OUTPUT_DIM = 155
    HID_DIM = 512
    ENC_LAYERS, DEC_LAYERS = 1, 1
    ENC_HEADS, DEC_HEADS = 32, 32
    ENC_PF_DIM, DEC_PF_DIM = HID_DIM*4, HID_DIM*4
    ENC_DROPOUT, DEC_DROPOUT = 0.1, 0.1

    enc = Transformer.Encoder(INPUT_DIM, HID_DIM, ENC_LAYERS, ENC_HEADS, ENC_PF_DIM, ENC_DROPOUT, device)
    dec = Transformer.Decoder(OUTPUT_DIM, HID_DIM, DEC_LAYERS, DEC_HEADS, DEC_PF_DIM, DEC_DROPOUT, device)

    tokenised_seq_PAD_IDX = 0
    TRG_PAD_IDX = 0

    model = Transformer.Seq2Seq(enc, dec, tokenised_seq_PAD_IDX, TRG_PAD_IDX, device).to(device)
    #mod_path = './model/shark_model.pt'
    mod_path = 'PLAbDab_nano/database_generate/shark_anarcii/model/shark_model.pt'
    model.load_state_dict(torch.load(mod_path, map_location=device))
    model.eval()

    # remove all B's and replace with X
    sequences = [seq.replace('B', 'X') for seq in ls]
    token_ls = tokenise_sharks(sequences) # Make the tokenised list of shark sequences for the model

    # condition at the following residues ============================================================
    sixty_seven = torch.tensor(Tokenizer.number_tokeniser(['67'])).unsqueeze(0)
    seventy_four = torch.tensor(Tokenizer.number_tokeniser(['74'])).unsqueeze(0)
    eighty_five = torch.tensor(Tokenizer.number_tokeniser(['85'])).unsqueeze(0)
    start = torch.tensor(Tokenizer.number_tokeniser(['2'])).unsqueeze(0)

    # Regex patterns we are looking for ==============================================================
    regex50_52 = re.compile(r'([NHPSK]E[ADQGEKV]|D[WE][EE])|KKE|NQE|NWE|NTE')
    regex99 = re.compile(r'SVNKGA')

    num_tokens = np.array(Tokenizer.num_vocab())
    aa_tokens = np.array(Tokenizer.seq_vocab())

    all_tuples = []
    with torch.no_grad():
        for tokenised_seq, raw_seq in zip(token_ls, sequences):
            batch_size = 1
            trg_len = tokenised_seq.shape[1]

            src_mask = model.make_src_mask(tokenised_seq)
            enc_src = model.encoder(tokenised_seq, src_mask)

            max_input = torch.zeros(batch_size, trg_len, device=device, dtype=torch.long)
            max_input[0, 0] = torch.tensor(1)
            input = torch.tensor(1).unsqueeze(0).unsqueeze(0) # start token of 1

            for t in range(1, trg_len):
                trg_mask = model.make_trg_mask(input)
                output, _ = model.decoder(input, enc_src, trg_mask, src_mask)

                # REMEMBER FOR THIS THAT FULL SEQ BEGINS AT 0
                if input[0][-1].item() in (51, 52, 53, 54, 55, 56) and regex50_52.search(raw_seq[t-3:t]):
                    pred_token = sixty_seven # modify the next token to be 68
                elif raw_seq[:3] == "ARV" and t == 1:
                    pred_token = start # modify the start token to be 2
                elif raw_seq[t-2:t] == "GG" and input[0][-1].item() == 75:
                    # print(raw_seq[t-2:t], input[0][-1].item())
                    pred_token = seventy_four # modify the start token to be 74 not 73
                elif regex99.search(raw_seq[t-7:t-1]):
                    print(input[0][-1].item(), raw_seq[t-7:t-1])
                    pred_token = eighty_five # modify the next token to be 85
                else:            
                    # Take last predicted output
                    pred_token = output.argmax(2)[:,-1].unsqueeze(1)
                
                max_input[:, t:t+1] = pred_token
                input = max_input[:, :t+1]
            
            max_input_tokens = num_tokens[max_input[:, :trg_len].to("cpu")]
            src_tokens = aa_tokens[tokenised_seq[:, :trg_len].to("cpu")]

            seq, num = [], []
            in_x_run, x_count = False, 0
            seq, num = [], []
            in_x_run, x_count = False, 0

            for seq_position in range(1,trg_len):
                if max_input_tokens[0, seq_position] == '<EOS>':
                    break
                elif max_input_tokens[0, seq_position] == 'X':
                    x_count += 1
                    in_x_run = True
                elif max_input_tokens[0, seq_position].isdigit() and in_x_run:
                    construction = build_inward_list(length=x_count, 
                                                    start_num=int(max_input_tokens[0, (seq_position - (x_count+1))]), # number before X began
                                                    end_num=int(max_input_tokens[0, seq_position])) # current number
                    
                    num[(seq_position - x_count):seq_position] = construction # add the construction over the previous sequence
                    num.append((max_input_tokens[0, seq_position], ' ')) # add the end
                    in_x_run=False
                    x_count = 0
                else:
                    num.append((max_input_tokens[0, seq_position], ' '))
                seq.append(src_tokens[0, seq_position])
            all_tuples.append(list(zip(num, seq)))
    
    return(all_tuples)






