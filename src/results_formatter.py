"""
Formatação e salvamento de resultados de execuções de algoritmos CSP.

Classes:
    ResultsFormatter: Armazena, formata e salva resultados detalhados e comparativos.

Métodos:
    add_algorithm_results(...): Adiciona resultados de um algoritmo.
    format_detailed_results(): Formata resultados detalhados.
    save_detailed_report(...): Salva relatório em arquivo.
"""
import statistics
from typing import List, Dict, Any
from tabulate import tabulate
from pathlib import Path
import logging

class ResultsFormatter:
    """
    Armazena, formata e salva resultados detalhados e comparativos de execuções de algoritmos CSP.

    Attributes:
        results (dict): Resultados por algoritmo.
        extra_info (dict): Informações extras do experimento.
    """
    def __init__(self):
        self.results = {}
        self.extra_info = {}
    
    def add_algorithm_results(self, algorithm_name: str, executions: List[Dict[str, Any]]):
        """
        Adiciona resultados de um algoritmo.

        Args:
            algorithm_name (str): Nome do algoritmo.
            executions (list): Lista de execuções detalhadas.
        """
        self.results[algorithm_name] = executions
    
    def format_detailed_results(self) -> str:
        """
        Formata os resultados detalhados para apresentação.

        Returns:
            str: Resultados formatados.
        """
        logger = logging.getLogger(__name__)
        output = []
        output.append("=" * 80)
        output.append("RELATÓRIO DETALHADO DE RESULTADOS")
        output.append("=" * 80)
        

        # Adicionar informações do dataset
        output.append(self._format_dataset_info())
        
        # Tabelas individuais por algoritmo
        for algorithm_name in self.results:
            output.append(self._format_algorithm_table(algorithm_name))
            output.append(self._format_algorithm_statistics(algorithm_name))
            output.append(self._format_algorithm_strings(algorithm_name))
            output.append("\n" + "=" * 80 + "\n")
        
        # Tabela comparativa final
        output.append(self._format_comparative_table())
        
        return "\n".join(output)

    def _format_algorithm_table(self, algorithm_name: str) -> str:
        """Formata tabela individual de um algoritmo"""
        executions = self.results[algorithm_name]
        output = [f"\n📊 TABELA DE EXECUÇÕES - {algorithm_name.upper()}"]
        output.append("-" * 60)
        
        # Cabeçalhos da tabela
        headers = [
            "Execução",
            "Tempo (s)",
            "Distância",
            "Status"
        ]
        # Dados da tabela
        table_data = []
        for i, exec_data in enumerate(executions, 1):
            distancia = exec_data.get('distancia', exec_data.get('melhor_distancia', '-'))
            status = "✓ OK"
            if exec_data.get('erro'):
                status = f"✗ {exec_data['erro']}"
                distancia = "-"
            elif exec_data.get('timeout'):
                status = "⏰ Timeout"
            elif distancia == float('inf'):
                status = "∞ Sem solução"
                distancia = "∞"
            row = [
                i,
                f"{exec_data['tempo']:.4f}",
                distancia,
                status
            ]
            table_data.append(row)
        table = tabulate(table_data, headers=headers, tablefmt="grid", stralign="center")
        output.append(table)
        return "\n".join(output)
    
    def _format_algorithm_statistics(self, algorithm_name: str) -> str:
        """Formata estatísticas detalhadas de um algoritmo"""
        executions = self.results[algorithm_name]
        valid_executions = [
            exec_data for exec_data in executions 
            if not exec_data.get('erro') and 
               exec_data.get('distancia', exec_data.get('melhor_distancia')) not in ['-', float('inf')]
        ]
        output = [f"\n📈 ESTATÍSTICAS DETALHADAS - {algorithm_name.upper()}"]
        output.append("-" * 60)
        
        if not valid_executions:
            output.append("❌ Nenhuma execução válida para calcular estatísticas.")
            return "\n".join(output)

        tempos = [exec_data['tempo'] for exec_data in valid_executions]
        distancias = [exec_data.get('distancia', exec_data.get('melhor_distancia', float('inf'))) 
                     for exec_data in valid_executions]

        stats_data = [
            ["TEMPO DE EXECUÇÃO", ""],
            ["Execuções Válidas", f"{len(valid_executions)}/{len(executions)}"],
            ["Média", f"{statistics.mean(tempos):.4f} s"],
            ["Mediana", f"{statistics.median(tempos):.4f} s"],
            ["Desvio Padrão", f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f} s"],
            ["Mínimo", f"{min(tempos):.4f} s"],
            ["Máximo", f"{max(tempos):.4f} s"],
        ]
        if distancias:
            stats_data.extend([
                ["", ""],
                ["DISTÂNCIA (em relação à string centro)", ""],
                ["Média", f"{statistics.mean(distancias):.2f}"],
                ["Mediana", f"{statistics.median(distancias):.2f}"],
                ["Desvio Padrão", f"{statistics.stdev(distancias) if len(distancias) > 1 else 0:.2f}"],
                ["Melhor (Mínima)", f"{min(distancias)}"],
                ["Pior (Máxima)", f"{max(distancias)}"],
            ])
        
        table = tabulate(stats_data, headers=["Métrica", "Valor"], tablefmt="grid", stralign="left")
        output.append(table)
        return "\n".join(output)
    
    def _format_algorithm_strings(self, algorithm_name: str) -> str:
        """Formata strings encontradas para auditoria"""
        executions = self.results[algorithm_name]
        output = [f"\n🔍 STRINGS PARA AUDITORIA - {algorithm_name.upper()}"]
        output.append("-" * 60)
        for i, exec_data in enumerate(executions, 1):
            distancia = exec_data.get('distancia', exec_data.get('melhor_distancia', '-'))
            string_result = exec_data.get('melhor_string', exec_data.get('string_resultado', ''))
            iteracoes = exec_data.get('iteracoes', exec_data.get('num_iteracoes', 0))
            output.append(f"Execução {i:2d}:")
            params = exec_data.get('params')
            seed = exec_data.get('seed')
            if params is not None:
                output.append(f"  Parâmetros de geração: {params}")
            if seed is not None:
                output.append(f"  Semente utilizada: {seed}")
            if exec_data.get('erro'):
                output.append(f"  ❌ Erro: {exec_data['erro']}")
                output.append(f"  Tempo: {exec_data['tempo']:.4f}s")
            else:
                output.append(f"  String Centro (Base): '{string_result}'")
                output.append(f"  Distância (algoritmo): {distancia}")
                output.append(f"  Iterações: {iteracoes}")
                output.append(f"  Tempo: {exec_data['tempo']:.4f}s")
            output.append("")
        return "\n".join(output)
    
    def _format_comparative_table(self) -> str:
        """Formata tabela comparativa entre algoritmos"""
        output = ["\n🏆 TABELA COMPARATIVA FINAL"]
        output.append("=" * 80)
        if not self.results:
            output.append("Nenhum resultado disponível para comparação.")
            return "\n".join(output)
        
        # Extrair dados comparativos
        comparative_data = []
        base_algorithms = {}
        
        for algorithm_full_name, executions in self.results.items():
            distancias = []
            tempos = []
            taxa_sucesso = 0.0
            
            parts = algorithm_full_name.split('_')
            if len(parts) >= 3 and 'Base' in parts[-2]:
                base_name = '_'.join(parts[:-2])
                base_num = parts[-2]
                algo_name = parts[-1]
                base_display = f"{base_name} ({base_num})"
            else:
                base_display = "N/A"
                algo_name = algorithm_full_name
                
            valid_executions = [
                exec_data for exec_data in executions 
                if not exec_data.get('erro') and 
                   exec_data.get('distancia', exec_data.get('melhor_distancia')) not in ['-', float('inf')]
            ]
            
            if not valid_executions:
                tempos = [exec_data['tempo'] for exec_data in executions]
                row = [
                    base_display,
                    algo_name,
                    f"{statistics.mean(tempos) if tempos else 0:.4f}",
                    f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}",
                    "ERRO",
                    "ERRO", 
                    "ERRO",
                    "0.0"
                ]
            else:
                tempos = [exec_data['tempo'] for exec_data in valid_executions]
                distancias = [exec_data.get('distancia', exec_data.get('melhor_distancia', float('inf'))) 
                             for exec_data in valid_executions]
                taxa_sucesso = len(valid_executions) / len(executions) * 100 if len(executions) > 0 else 0
                
                row = [
                    base_display,
                    algo_name,
                    f"{statistics.mean(tempos) if tempos else 0:.4f}",
                    f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}",
                    f"{min(distancias) if distancias else 'ERRO'}",
                    f"{statistics.mean(distancias):.2f}" if distancias else "ERRO",
                    f"{statistics.stdev(distancias) if len(distancias) > 1 else 0:.2f}" if distancias else "ERRO",
                    f"{taxa_sucesso:.1f}"
                ]
                
            comparative_data.append(row)
            
            if base_display not in base_algorithms:
                base_algorithms[base_display] = []
                
            if valid_executions:
                base_algorithms[base_display].append({
                    'name': algo_name,
                    'dist': min(distancias) if distancias else float('inf'),
                    'time': statistics.mean(tempos) if tempos else 0,
                    'success_rate': taxa_sucesso if isinstance(taxa_sucesso, (int, float)) else 0
                })
            else:
                base_algorithms[base_display].append({
                    'name': algo_name,
                    'dist': float('inf'),
                    'time': statistics.mean(tempos) if tempos else 0,
                    'success_rate': 0
                })
        
        headers = [
            "Base",
            "Algoritmo",
            "Tempo Médio (s)",
            "Desvio Tempo",
            "Melhor Distância",
            "Distância Média", 
            "Desvio Distancia",
            "Taxa Sucesso (%)"
        ]
        
        def sort_key(row):
            base_name = row[0]
            if row[4] == "ERRO":
                dist = float('inf')
            else:
                try:
                    dist = float(row[4])
                except Exception:
                    dist = float('inf')
            return (base_name, dist, float(row[2]))
        
        comparative_data.sort(key=sort_key)
        table = tabulate(comparative_data, headers=headers, tablefmt="grid", stralign="center")
        output.append(table)
        
        # Adicionar ranking POR BASE
        output.append("\n🥇 RANKING POR BASE:")
        output.append("-" * 60)
        
        for base_name, algorithms in sorted(base_algorithms.items()):
            output.append(f"\n💠 BASE: {base_name}")
            output.append("-" * 40)
            
            sorted_algorithms = sorted(algorithms, key=lambda x: (x['dist'], x['time']))
            
            for i, alg in enumerate(sorted_algorithms, 1):
                if alg['dist'] == float('inf'):
                    medal = "❌"
                    info = f"Falha na execução | Tempo: {alg['time']:.4f}s"
                else:
                    medal = "🥇" if i == 1 else "🥈" if i == 2 else "🥉" if i == 3 else f"{i}°"
                    info = f"Distância: {alg['dist']} | Tempo: {alg['time']:.4f}s | Sucesso: {alg['success_rate']:.1f}%"
                output.append(f"{medal} {alg['name']} - {info}")
                
        return "\n".join(output)
    
    def save_detailed_report(self, filename: str):
        """
        Salva o relatório detalhado em um arquivo.
        
        Args:
            filename: Nome do arquivo ou caminho completo
        """
        file_path = Path(filename)
        # Garantir que o diretório existe
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(self.format_detailed_results())
    
    def get_detailed_report_data(self):
        """
        Obtém os dados do relatório detalhado.

        Returns:
            dict: Dados do relatório, incluindo informações extras e resultados.
        """
        return {
            'extra_info': self.extra_info,
            'results': self.results
        }
    
    def _format_individual_report(self, algorithm_name, executions, base_name=None):
        """Formata relatório individual para um algoritmo específico"""
        output = []
        
        # Cabeçalho com destaque para algoritmo
        output.append(f"\n📊 RESULTADOS PARA {algorithm_name}")
        output.append("=" * 80)
        
        # Adicionar informações sobre a base/dataset
        if base_name:
            output.append(f"Base de Dados: {base_name}")
            
        # Verificar se há informações específicas sobre o centro/base string
        if hasattr(self, 'extra_info') and self.extra_info:
            # Outras informações opcionais
            seed = self.extra_info.get('seed')
            if seed:
                output.append(f"Seed: {seed}")
                
    def _format_dataset_info(self) -> str:
        """Formata informações do dataset incluindo string base de geração"""
        output = ["\n📊 INFORMAÇÕES DO DATASET"]
        output.append("-" * 60)
        
        if not hasattr(self, 'extra_info') or not self.extra_info:
            output.append("❌ Nenhuma informação do dataset disponível.")
            return "\n".join(output)
        
        # Extrair informações do dataset
        dataset_info = {}
        base_strings = set()
        
        for key, info in self.extra_info.items():
            if isinstance(info, dict):
                params = info.get('params', {})
                if params:
                    # Atualizar informações do dataset
                    for param_key, param_value in params.items():
                        if param_key not in dataset_info:
                            dataset_info[param_key] = param_value
                    
                    # Coletar string base se disponível
                    base_string = params.get('base_string')
                    if base_string:
                        base_strings.add(base_string)
        
        # Exibir informações básicas do dataset
        if dataset_info:
            output.append("📋 Parâmetros do Dataset:")
            for key, value in dataset_info.items():
                if key != 'base_string':  # Trataremos a base_string separadamente
                    if key == 'n':
                        output.append(f"  • Número de strings: {value}")
                    elif key == 'L':
                        output.append(f"  • Comprimento das strings: {value}")
                    elif key == 'alphabet':
                        output.append(f"  • Alfabeto: '{value}' (|Σ| = {len(str(value))})")
                    elif key == 'noise':
                        if isinstance(value, list):
                            output.append(f"  • Taxa de ruído: variável (min: {min(value):.3f}, max: {max(value):.3f})")
                        else:
                            output.append(f"  • Taxa de ruído: {value}")
                    elif key == 'seed':
                        output.append(f"  • Semente: {value}")
                    elif key == 'fully_random':
                        output.append(f"  • Modo: {'Totalmente aleatório' if value else 'Base + ruído'}")
        
        # Exibir string(s) base para auditoria
        if base_strings:
            output.append("")
            output.append("🎯 STRING BASE PARA AUDITORIA:")
            for i, base_string in enumerate(sorted(base_strings), 1):
                output.append(f"  Base {i}: ")
                output.append(f"    • String: '{base_string}'")
                dist = self.extra_info.get('distancia_string_base', '')
                output.append(f"    • Distância da string base: {dist}")
        else:
            output.append("❌ String base não disponível (dataset não sintético ou informação não capturada).")
        
        return "\n".join(output)