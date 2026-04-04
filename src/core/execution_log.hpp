/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/execution_log.hpp                                        */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Espelha std::cout e std::cerr para um arquivo de log por caso.  */
/*****************************************************************************/
/* Observacao: A ideia e preservar a experiencia normal no terminal enquanto  */
/* a mesma trilha textual fica registrada em run.log para auditoria didatica. */
/*****************************************************************************/

#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <streambuf>
#include <string>

namespace execution_log
{

/******************************************************************************/
/* CLASSE: TeeStreamBuf                                                       */
/* DESCRICAO: Stream buffer simples que duplica toda escrita para dois        */
/* destinos: o fluxo original do terminal e o arquivo de log do caso.        */
/******************************************************************************/
class TeeStreamBuf : public std::streambuf
{
public:
    TeeStreamBuf() = default;

    void set_primary(std::streambuf *primary)
    {
        primary_ = primary;
    }

    void set_secondary(std::streambuf *secondary)
    {
        secondary_ = secondary;
    }

protected:
    int overflow(int ch) override
    {
        if (traits_type::eq_int_type(ch, traits_type::eof()))
        {
            return sync() == 0 ? traits_type::not_eof(ch) : traits_type::eof();
        }

        const char c = traits_type::to_char_type(ch);
        const int primary_result =
            primary_ ? primary_->sputc(c) : traits_type::not_eof(ch);
        const int secondary_result =
            secondary_ ? secondary_->sputc(c) : traits_type::not_eof(ch);

        if ((primary_ && traits_type::eq_int_type(primary_result, traits_type::eof())) ||
            (secondary_ && traits_type::eq_int_type(secondary_result, traits_type::eof())))
        {
            return traits_type::eof();
        }

        return traits_type::not_eof(ch);
    }

    std::streamsize xsputn(const char *s, std::streamsize count) override
    {
        const std::streamsize primary_written =
            primary_ ? primary_->sputn(s, count) : count;
        const std::streamsize secondary_written =
            secondary_ ? secondary_->sputn(s, count) : count;
        return std::min(primary_written, secondary_written);
    }

    int sync() override
    {
        const int primary_sync = primary_ ? primary_->pubsync() : 0;
        const int secondary_sync = secondary_ ? secondary_->pubsync() : 0;
        return (primary_sync == 0 && secondary_sync == 0) ? 0 : -1;
    }

private:
    std::streambuf *primary_ = nullptr;
    std::streambuf *secondary_ = nullptr;
};

/******************************************************************************/
/* CLASSE: ExecutionLogScope                                                  */
/* DESCRICAO: RAII que instala um tee em cout/cerr para espelhar toda a       */
/* execucao textual em run.log, restaurando os buffers originais ao final.    */
/******************************************************************************/
class ExecutionLogScope
{
public:
    explicit ExecutionLogScope(const std::string &log_path)
        : file_path_(log_path)
    {
        old_cout_ = std::cout.rdbuf();
        old_cerr_ = std::cerr.rdbuf();

        log_file_.open(file_path_, std::ios::out | std::ios::trunc);
        if (!log_file_)
        {
            error_message_ = "falha ao abrir arquivo para escrita";
            return;
        }

        cout_tee_.set_primary(old_cout_);
        cout_tee_.set_secondary(log_file_.rdbuf());
        cerr_tee_.set_primary(old_cerr_);
        cerr_tee_.set_secondary(log_file_.rdbuf());

        std::cout.rdbuf(&cout_tee_);
        std::cerr.rdbuf(&cerr_tee_);
        active_ = true;
    }

    ExecutionLogScope(const ExecutionLogScope &) = delete;
    ExecutionLogScope &operator=(const ExecutionLogScope &) = delete;

    ~ExecutionLogScope()
    {
        if (!active_)
            return;

        std::cout.flush();
        std::cerr.flush();
        std::cout.rdbuf(old_cout_);
        std::cerr.rdbuf(old_cerr_);
        log_file_.flush();
    }

    bool active() const
    {
        return active_;
    }

    const std::string &file_path() const
    {
        return file_path_;
    }

    const std::string &error_message() const
    {
        return error_message_;
    }

private:
    std::string file_path_;
    std::string error_message_;
    std::ofstream log_file_;
    TeeStreamBuf cout_tee_;
    TeeStreamBuf cerr_tee_;
    std::streambuf *old_cout_ = nullptr;
    std::streambuf *old_cerr_ = nullptr;
    bool active_ = false;
};

} // namespace execution_log
